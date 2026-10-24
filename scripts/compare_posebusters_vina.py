"""Compare a Vina run with the PoseBusters authors' own Vina results, complex by complex.

The PoseBusters Zenodo record (8278563) publishes per-complex results for every method in
`posebusters_paper_results.csv`. This script downloads that file (once, into data/raw/) and
compares its Vina top-1 outcomes, without post-processing, with one of our runs.

Usage:
    uv run python scripts/compare_posebusters_vina.py [RUN_DIR]

RUN_DIR defaults to results/vina/baseline.
"""

import sys
import urllib.request
from pathlib import Path

import pandas as pd

from dockprepbench.analysis.stats import compare, load_run, success_rate, successes

URL = "https://zenodo.org/api/records/8278563/files/posebusters_paper_results.csv/content"
LOCAL = Path("data/raw/posebusters_paper_results.csv")
CHECKS_PREFIXES = ("docked_", "molecule_", "molecular_", "sp3_", "double_", "bond_", "no_",
                   "aromatic_", "energy_")  # fmt: skip


def their_vina() -> pd.DataFrame:
    """
    Load the PoseBusters authors' Vina results on the PoseBusters Benchmark set.

    Returns:
        pd.DataFrame: Indexed by complex id (PDB_CCD), with `top1` and `pb_valid` flags.
    """
    if not LOCAL.is_file():
        LOCAL.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(URL, LOCAL)
    pb = pd.read_csv(LOCAL)
    vina = pb[(pb.dataset == "posebuster") & (pb.method == "vina")]
    vina = vina[vina["post-processing"] == "none"]
    checks = [c for c in pb.columns if c.startswith(CHECKS_PREFIXES)]
    return pd.DataFrame(
        {
            "top1": (vina["rmsd"] < 2).to_numpy(),
            "pb_valid": (
                (vina["rmsd"] < 2) & vina[checks].fillna(False).astype(bool).all(axis=1)
            ).to_numpy(),
        },
        index=vina["pdb_id"] + "_" + vina["ccd_id"],
    )


def main() -> None:
    """Print rates, the paired comparison, and the same on the complexes our run docked."""
    run_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results/vina/baseline")
    run = load_run(run_dir / "results.csv")
    theirs = their_vina()
    ours = successes(run)
    shared = ours.index.intersection(theirs.index)
    print(f"PoseBusters Vina: top-1 {success_rate(theirs.loc[shared, 'top1'])}")
    print(f"                  PB-valid {success_rate(theirs.loc[shared, 'pb_valid'])}")
    print(f"{run_dir}: top-1 {success_rate(ours[shared])}")
    for label, ids in [
        ("all", shared),
        ("docked by us", shared.intersection(run.index[run.error == ""])),
    ]:
        c = compare(theirs.loc[ids, "top1"], ours[ids])
        print(
            f"{label} ({len(ids)}): theirs {100 * theirs.loc[ids, 'top1'].mean():.1f}%, "
            f"ours {100 * ours[ids].mean():.1f}%; only theirs {len(c.only_a)}, only ours "
            f"{len(c.only_b)}; exact McNemar p = {c.p_value:.3g}"
        )


if __name__ == "__main__":
    main()
