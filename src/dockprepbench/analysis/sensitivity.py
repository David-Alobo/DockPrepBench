"""Which complexes change outcome with a preparation choice: flips against ligand properties.

For each engine and arm, the complexes whose top-1 outcome differs from the baseline
("flipped") are compared with those whose outcome does not, on properties of the ligand
and its complex. The comparisons are exploratory: they describe where flips concentrate,
not why.
"""

from collections.abc import Iterable

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from scipy.stats import fisher_exact, mannwhitneyu

from dockprepbench.analysis.report import EngineRuns
from dockprepbench.analysis.stats import holm_adjust, successes
from dockprepbench.datasets.posebusters import Complex
from dockprepbench.io import load_sdf_file
from dockprepbench.preparation.ligand import protonate_ligand

# Crystallisation additives, solvent and simple ions: present in many protein files but not
# cofactors. Residue names as in the PDB chemical component dictionary.
ADDITIVES = frozenset({
    "HOH", "EDO", "GOL", "SO4", "PO4", "DMS", "PEG", "PG4", "ACT", "CL", "NA", "K", "BR",
    "IOD", "FMT", "MPD", "TRS", "EPE", "MES", "BME", "IMD", "NO3", "SCN", "1PE",
})  # fmt: skip
NUMERIC = ("heavy_atoms", "rotatable_bonds", "abs_charge")
BINARY = ("charge_changed", "cofactor", "seed_sensitive")


def cofactors(entry: Complex) -> list[str]:
    """
    List the cofactors in a complex's protein file.

    Args:
        entry (Complex): Benchmark complex.

    Returns:
        list[str]: Residue names of heterogens other than additives and the ligand of
            interest (whose CCD code ends the complex ID), sorted.
    """
    hets = {
        line[17:20].strip()
        for line in entry.protein.read_text().splitlines()
        if line.startswith("HETATM")
    }
    return sorted(hets - ADDITIVES - {entry.id.split("_")[-1]})


def ligand_properties(complexes: Iterable[Complex], ph: float = 7.4) -> pd.DataFrame:
    """
    Describe each complex's ligand as docked in the baseline, and its protein file.

    Args:
        complexes (Iterable[Complex]): Benchmark complexes.
        ph (float): pH at which protonation is assigned in the protonation arm.

    Returns:
        pd.DataFrame: Indexed by complex id: heavy_atoms, rotatable_bonds, charge and
            abs_charge of the supplied ligand, charge_changed (protonation at `ph` changes
            the formal charge) and cofactor (any cofactor in the protein file).
    """
    rows = []
    for entry in complexes:
        ligand = load_sdf_file(entry.start_ligand)[0]
        charge = Chem.GetFormalCharge(ligand)
        try:
            changed = Chem.GetFormalCharge(protonate_ligand(ligand, ph)) != charge
        except Exception:  # noqa: BLE001 - a ligand Dimorphite cannot handle is "unknown"
            changed = None
        rows.append({
            "id": entry.id,
            "heavy_atoms": ligand.GetNumHeavyAtoms(),
            "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(ligand),
            "charge": charge,
            "abs_charge": abs(charge),
            "charge_changed": changed,
            "cofactor": bool(cofactors(entry)),
        })  # fmt: skip
    return pd.DataFrame(rows).set_index("id")


def seed_sensitive(runs: EngineRuns) -> pd.Series | None:
    """
    Flag complexes whose top-1 outcome differs between seed replicates of the baseline.

    Args:
        runs (EngineRuns): One engine's runs.

    Returns:
        pd.Series | None: Boolean per complex, or None without seed replicates.
    """
    if not runs.seeds:
        return None
    table = pd.DataFrame(
        {k: successes(v) for k, v in {"baseline": runs.baseline, **runs.seeds}.items()}
    )
    return ~(table.all(axis=1) | (~table).all(axis=1))


def _compare(flipped: pd.Series, values: pd.Series, binary: bool) -> dict:
    """
    Compare a property between flipped and unflipped complexes.

    Args:
        flipped (pd.Series): Boolean per complex.
        values (pd.Series): Property per complex, aligned with `flipped`; NaN is dropped.
        binary (bool): Fisher's exact test on a boolean property; otherwise Mann-Whitney U.

    Returns:
        dict: Group summaries, an effect size and the two-sided p-value.
    """
    keep = values.notna()
    flipped, values = flipped[keep], values[keep]
    yes, no = values[flipped], values[~flipped]
    if binary:
        values = values.astype(bool)
        table = [
            [int((values & flipped).sum()), int((values & ~flipped).sum())],
            [int((~values & flipped).sum()), int((~values & ~flipped).sum())],
        ]
        odds, p = fisher_exact(table)
        with_flip = table[0][0] / max(1, table[0][0] + table[0][1])
        without_flip = table[1][0] / max(1, table[1][0] + table[1][1])
        return {
            "flipped_value": float(yes.astype(bool).mean()),
            "unflipped_value": float(no.astype(bool).mean()),
            "summary": "share with property",
            "effect": float(odds),
            "effect_kind": "odds ratio",
            "flip_rate_with": with_flip,
            "flip_rate_without": without_flip,
            "p": float(p),
        }
    u, p = mannwhitneyu(yes, no, alternative="two-sided")
    return {
        "flipped_value": float(yes.median()),
        "unflipped_value": float(no.median()),
        "summary": "median",
        # Rank-biserial correlation: +1 when every flipped complex has the larger value.
        "effect": float(2 * u / (len(yes) * len(no)) - 1),
        "effect_kind": "rank-biserial r",
        "flip_rate_with": float("nan"),
        "flip_rate_without": float("nan"),
        "p": float(p),
    }


def sensitivity_table(runs: dict[str, EngineRuns], properties: pd.DataFrame) -> pd.DataFrame:
    """
    Compare flipped with unflipped complexes on each property, per engine and arm.

    Args:
        runs (dict[str, EngineRuns]): Output of `report.discover`.
        properties (pd.DataFrame): Output of `ligand_properties`.

    Returns:
        pd.DataFrame: One row per engine, arm and property, with p-values Holm-adjusted
            over the properties tested for that engine and arm.
    """
    rows = []
    for e in runs.values():
        props = properties.copy()
        sensitive = seed_sensitive(e)
        if sensitive is not None:
            props["seed_sensitive"] = sensitive
        base = successes(e.baseline)
        for arm, run in e.arms.items():
            hits = successes(run)
            shared = base.index.intersection(hits.index).intersection(props.index)
            flipped = base[shared] != hits[shared]
            if flipped.all() or not flipped.any():
                continue
            family = []
            for prop in (*NUMERIC, *BINARY):
                if prop not in props:
                    continue
                result = _compare(flipped, props.loc[shared, prop], binary=prop in BINARY)
                family.append({
                    "engine": e.engine,
                    "arm": arm,
                    "property": prop,
                    "n_flipped": int(flipped.sum()),
                    "n": len(shared),
                } | result)  # fmt: skip
            for row, p_holm in zip(family, holm_adjust([r["p"] for r in family]), strict=True):
                rows.append(row | {"p_holm": p_holm})
    return pd.DataFrame(rows)
