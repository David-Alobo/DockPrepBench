"""Receptor preparation: repair and protonation."""

import io
from pathlib import Path

from openmm import Context, LocalEnergyMinimizer, Platform, VerletIntegrator, unit
from openmm.app import CutoffNonPeriodic, ForceField, Modeller, PDBFile, element
from openmm.app.topology import Atom
from pdbfixer import PDBFixer

DEFAULT_PH = 7.4


def protonate_receptor(pdb_path: str | Path, ph: float = DEFAULT_PH) -> str:
    """
    Prepare a protein for docking with PDBFixer.

    Nonstandard residues are mutated to their standard parents, missing heavy atoms are
    rebuilt and hydrogens are added for the given pH. The added atoms are then relaxed
    with the crystal atoms held fixed. Missing loops are not modelled, and all heterogens,
    including waters and cofactors, are removed.

    Args:
        pdb_path (str | Path): Input protein PDB file.
        ph (float): pH used to assign protonation states.

    Returns:
        str: Prepared protein as PDB text, keeping the original chain and residue IDs.
    """
    fixer = PDBFixer(filename=str(pdb_path))
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    _remove_hydrogens(fixer)
    crystal_atoms = {_atom_key(atom) for atom in fixer.topology.atoms()}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=ph)
    positions = _relax_added_atoms(fixer, crystal_atoms)

    out = io.StringIO()
    PDBFile.writeFile(fixer.topology, positions, out, keepIds=True)
    return out.getvalue()


def _remove_hydrogens(fixer: PDBFixer) -> None:
    """
    Delete any hydrogens already in the input, so protonation follows this protocol only.

    Input hydrogens may also carry naming conventions the force field templates reject.

    Args:
        fixer (PDBFixer): Fixer to modify in place.
    """
    modeller = Modeller(fixer.topology, fixer.positions)
    modeller.delete([a for a in modeller.topology.atoms() if a.element == element.hydrogen])
    fixer.topology, fixer.positions = modeller.topology, modeller.positions


def _relax_added_atoms(fixer: PDBFixer, fixed: set[tuple]) -> unit.Quantity:
    """
    Energy-minimise atoms added by PDBFixer while keeping the given atoms in place.

    PDBFixer places rebuilt side chains, terminal oxygens and hydrogens crudely, sometimes
    close enough to neighbours that distance-based bond perception (as in meeko) finds
    bonds that do not exist. A short Amber14 minimisation fixes the geometry.

    Args:
        fixer (PDBFixer): Fixer after atoms and hydrogens have been added.
        fixed (set[tuple]): Keys (see `_atom_key`) of atoms that must not move.

    Returns:
        unit.Quantity: Minimised positions; fixed atoms keep their input coordinates.
    """
    system = ForceField("amber14-all.xml").createSystem(
        fixer.topology, nonbondedMethod=CutoffNonPeriodic, nonbondedCutoff=1.0 * unit.nanometer
    )
    for atom in fixer.topology.atoms():
        if _atom_key(atom) in fixed:
            system.setParticleMass(atom.index, 0)  # OpenMM never moves massless particles

    context = Context(system, VerletIntegrator(0.001), Platform.getPlatformByName("CPU"))
    context.setPositions(fixer.positions)
    LocalEnergyMinimizer.minimize(context, tolerance=10.0, maxIterations=500)
    return context.getState(getPositions=True).getPositions()


def _atom_key(atom: Atom) -> tuple[str, str, str, str]:
    """
    Identify an atom independently of its index, which changes as PDBFixer adds atoms.

    Args:
        atom (Atom): OpenMM topology atom.

    Returns:
        tuple[str, str, str, str]: Chain ID, residue number, insertion code and atom name.
    """
    residue = atom.residue
    return (residue.chain.id, residue.id, residue.insertionCode, atom.name)
