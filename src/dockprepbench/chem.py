"""Small RDKit helpers shared across stages."""

from rdkit import Chem


def skeleton(mol: Chem.Mol) -> Chem.Mol:
    """
    Strip a molecule to heavy atoms with no formal charges, keeping bonds and coordinates.

    Used to match atoms between two forms of the same molecule that differ only in
    protonation state. Atom indices equal those of `Chem.RemoveAllHs(mol)`.

    Args:
        mol (Chem.Mol): Molecule, with or without a conformer.

    Returns:
        Chem.Mol: Charge-free heavy-atom copy, for atom matching only.
    """
    stripped = Chem.RWMol(Chem.RemoveAllHs(mol))
    for atom in stripped.GetAtoms():
        atom.SetFormalCharge(0)
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(True)
    # Not re-sanitised: neutral valences can be invalid (e.g. nitro N) and aromatic
    # rings can lose the H they need to kekulize; matching only needs elements, bond
    # orders and aromaticity flags, which are unchanged.
    stripped.UpdatePropertyCache(strict=False)
    return stripped.GetMol()
