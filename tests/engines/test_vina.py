import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from dockprepbench.engines.base import Box, DockingResult
from dockprepbench.engines.vina import VinaEngine
from dockprepbench.io import load_pdb_file

pytest.importorskip("vina", reason="vina bindings are only available on Linux x86_64")


def test_vina_docks_into_box(peptide_pdb):
    atoms = list(load_pdb_file(peptide_pdb).get_atoms())
    center = tuple(sum(a.coord[i] for a in atoms) / len(atoms) for i in range(3))
    box = Box(center=center, size=(20.0, 20.0, 20.0))
    ligand = Chem.AddHs(Chem.MolFromSmiles("NC(=O)c1ccccc1"))
    AllChem.EmbedMolecule(ligand, randomSeed=1)

    result = VinaEngine(exhaustiveness=1, n_poses=3).dock(peptide_pdb, ligand, box)

    assert isinstance(result, DockingResult)
    assert 1 <= len(result.poses) <= 3
    assert result.scores == sorted(result.scores)
    for pose in result.poses:
        assert Chem.MolToSmiles(Chem.RemoveHs(pose)) == "NC(=O)c1ccccc1"
        # Vina confines heavy atoms; meeko re-adds nonpolar hydrogens that may poke out.
        heavy = Chem.RemoveHs(pose).GetConformer().GetPositions()
        for x, c, s in zip(heavy.T, box.center, box.size, strict=True):
            assert (abs(x - c) <= s / 2).all()


def test_vina_is_reproducible_with_fixed_seed(peptide_pdb):
    box = Box(center=(0.0, 0.0, 0.0), size=(20.0, 20.0, 20.0))
    ligand = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(ligand, randomSeed=1)
    engine = VinaEngine(exhaustiveness=1, n_poses=1)

    first = engine.dock(peptide_pdb, ligand, box)
    second = engine.dock(peptide_pdb, ligand, box)

    assert first.scores == second.scores
