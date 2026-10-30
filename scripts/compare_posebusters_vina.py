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
    """
    Print rates and the paired comparison (all complexes, and those our run docked).

    With a second argument, also write the comparison as a CSV table for the paper.
    """
    run_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results/vina/baseline")
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else None
    run = load_run(run_dir / "results.csv")
    theirs = their_vina()
    ours = successes(run)
    shared = ours.index.intersection(theirs.index)
    print(f"PoseBusters Vina: top-1 {success_rate(theirs.loc[shared, 'top1'])}")
    print(f"                  PB-valid {success_rate(theirs.loc[shared, 'pb_valid'])}")
    print(f"{run_dir}: top-1 {success_rate(ours[shared])}")
    failed = shared.intersection(run.index[run.error != ""])
    rows = []
    for label, ids in [
        ("all", shared),
        ("docked", shared.intersection(run.index[run.error == ""])),
    ]:
        c = compare(theirs.loc[ids, "top1"], ours[ids])
        rate = success_rate(theirs.loc[ids, "top1"])
        print(
            f"{label} ({len(ids)}): theirs {100 * theirs.loc[ids, 'top1'].mean():.1f}%, "
            f"ours {100 * ours[ids].mean():.1f}%; only theirs {len(c.only_a)}, only ours "
            f"{len(c.only_b)}; exact McNemar p = {c.p_value:.3g}"
        )
        rows.append({
            "set": label,
            "n": len(ids),
            "theirs_top1": rate.value,
            "theirs_top1_low": rate.low,
            "theirs_top1_high": rate.high,
            "theirs_pb_valid": theirs.loc[ids, "pb_valid"].mean(),
            "ours_top1": ours[ids].mean(),
            "diff_pts": 100 * (theirs.loc[ids, "top1"].mean() - ours[ids].mean()),
            "only_theirs": len(c.only_a),
            "only_ours": len(c.only_b),
            "agree": c.n - len(c.only_a) - len(c.only_b),
            "p": c.p_value,
            "theirs_solved_our_failures": int(theirs.loc[failed, "top1"].sum()),
        })  # fmt: skip
    if out:
        pd.DataFrame(rows).to_csv(out, index=False, lineterminator="\n")
        print(f"Wrote {out}")


if __name__ == "__main__":
    main()
