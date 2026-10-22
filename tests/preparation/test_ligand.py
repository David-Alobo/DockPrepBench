import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from dockprepbench.evaluation.rmsd import pose_rmsd
from dockprepbench.preparation.ligand import generate_conformer, prepare_ligand, protonate_ligand


def embedded(smiles, seed=1):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    return mol


def charges(mol):
    return sorted(a.GetFormalCharge() for a in mol.GetAtoms() if a.GetFormalCharge())


def test_protonation_charges_acids_and_amines():
    # beta-alanine: neutral input, zwitterion at pH 7.4
    mol = embedded("NCCC(=O)O")

    protonated = protonate_ligand(mol, ph=7.4)

    assert charges(protonated) == [-1, 1]
    assert Chem.MolToSmiles(Chem.RemoveHs(protonated)) == "[NH3+]CCC(=O)[O-]"


def test_protonation_keeps_heavy_atom_coordinates():
    mol = embedded("NCCC(=O)O")

    protonated = protonate_ligand(mol, ph=7.4)

    assert pose_rmsd(protonated, mol) == pytest.approx(0.0, abs=1e-6)


def test_protonation_adds_hydrogens_with_coordinates():
    protonated = protonate_ligand(embedded("CCN"), ph=7.4)

    positions = protonated.GetConformer().GetPositions()
    assert protonated.GetNumAtoms() == Chem.AddHs(Chem.MolFromSmiles("CC[NH3+]")).GetNumAtoms()
    assert not np.isclose(positions, 0.0).all(axis=1).any()


def test_protonation_leaves_neutral_ligand_unchanged():
    mol = embedded("NC(=O)c1ccccc1")

    protonated = protonate_ligand(mol, ph=7.4)

    assert charges(protonated) == []
    assert Chem.MolToSmiles(Chem.RemoveHs(protonated)) == Chem.MolToSmiles(Chem.RemoveHs(mol))


def test_generate_conformer_is_seeded_and_new():
    mol = embedded("CCOc1ccccc1C(=O)N", seed=99)

    first = generate_conformer(mol, seed=0)
    again = generate_conformer(mol, seed=0)

    assert np.allclose(first.GetConformer().GetPositions(), again.GetConformer().GetPositions())
    assert pose_rmsd(first, mol) > 0.1
    assert first.GetNumConformers() == 1


def test_prepare_ligand_supplied_is_identity():
    mol = embedded("CCN")

    assert prepare_ligand(mol, "supplied", "supplied", ph=7.4) is mol


def test_prepare_ligand_protonates_before_new_conformer():
    prepared = prepare_ligand(embedded("CCC(=O)O"), "dimorphite", "rdkit", ph=7.4)

    assert charges(prepared) == [-1]
    assert prepared.GetNumConformers() == 1


def test_prepare_ligand_rejects_unknown_choice():
    with pytest.raises(ValueError, match="protonation"):
        prepare_ligand(embedded("CCN"), "obabel", "supplied", ph=7.4)
