"""Follow-up Vina arms: protonation at each box padding, and a box sized to the ligand.

These runs test assumptions of the one-factor-at-a-time design and are analysed apart from the
primary arms, each set with its own correction (paper section on the interaction).

1. Protonation x box: ligand protonation at pH 7.4 repeated at 2, 6 and 8 Ã… padding, compared
   with the supplied ligand at the same padding (4 Ã… is the primary protonation arm). Holm over
   the four paddings. Whether the effect differs from the effect at the baseline padding is
   tested complex by complex with a sign test on the difference in effect.
2. The box from the ligand's radius of gyration (Feinstein and Brylinski 2015) against the
   baseline and the 2 Ã… box, with box volumes and whether the crystal ligand fits inside.

Run with `python -m dockprepbench.analysis.followup RESULTS DATASET OUT`.
"""

import argparse
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from scipy.stats import binomtest

from dockprepbench.analysis.stats import (
    compare,
    holm_adjust,
    load_run,
    paired_difference_interval,
    success_rate,
    successes,
)
from dockprepbench.datasets.posebusters import load_complexes
from dockprepbench.io import load_sdf_file
from dockprepbench.preparation.box import DEFAULT_PADDING, box_from_gyration, box_from_ligand

# Padding in Ã… -> (supplied-ligand run, protonated run).
PADDINGS = {
    2: ("box_padding_2", "ligand_dimorphite_box_2"),
    4: ("baseline", "ligand_dimorphite"),
    6: ("box_padding_6", "ligand_dimorphite_box_6"),
    8: ("box_padding_8", "ligand_dimorphite_box_8"),
}


def _pair(a: pd.Series, b: pd.Series) -> dict:
    """
    Paired comparison of run b against run a: rates, difference with interval, counts, p.

    Args:
        a (pd.Series): Success flags of the reference run.
        b (pd.Series): Success flags of the compared run.

    Returns:
        dict: rate_a, low_a, high_a, rate_b, low_b, high_b, diff_pts, diff_low, diff_high,
            gained, lost, p.
    """
    c = compare(a, b)
    ra, rb = (
        success_rate(a.loc[a.index.intersection(b.index)]),
        success_rate(b.loc[a.index.intersection(b.index)]),
    )
    low, high = paired_difference_interval(c.both, len(c.only_a), len(c.only_b), c.n)
    return {
        "rate_a": ra.value, "low_a": ra.low, "high_a": ra.high,
        "rate_b": rb.value, "low_b": rb.low, "high_b": rb.high,
        "diff_pts": 100 * (rb.value - ra.value), "diff_low": 100 * low, "diff_high": 100 * high,
        "gained": len(c.only_b), "lost": len(c.only_a), "p": c.p_value,
    }  # fmt: skip


def interaction_table(vina: Path) -> pd.DataFrame:
    """
    Protonation against the supplied ligand at each box padding, Holm over the four paddings.

    Each row also tests, complex by complex, whether the protonation effect at that padding
    differs from the effect at the baseline padding (sign test on the difference in effect),
    Holm over the three non-baseline paddings.

    Args:
        vina (Path): Vina results directory holding the primary and follow-up runs.

    Returns:
        pd.DataFrame: One row per padding.
    """
    hits = {run: successes(load_run(vina / run / "results.csv"))
            for pair in PADDINGS.values() for run in pair}  # fmt: skip
    effect = {
        pad: hits[prot].astype(int) - hits[sup].astype(int) for pad, (sup, prot) in PADDINGS.items()
    }
    rows = []
    for pad, (sup, prot) in PADDINGS.items():
        row = {"padding": pad, **_pair(hits[sup], hits[prot])}
        if pad != 4:
            d = (effect[pad] - effect[4]).dropna()
            up, down = int((d > 0).sum()), int((d < 0).sum())
            row["vs_baseline_padding_p"] = (
                binomtest(up, up + down, 0.5).pvalue if up + down else 1.0
            )
        rows.append(row)
    table = pd.DataFrame(rows)
    table["p_holm"] = holm_adjust(table["p"].tolist())
    others = table["vs_baseline_padding_p"].notna()
    table.loc[others, "vs_baseline_padding_p_holm"] = holm_adjust(
        table.loc[others, "vs_baseline_padding_p"].tolist()
    )
    return table


def combined_table(vina: Path) -> pd.DataFrame:
    """
    The combined protocol (2 Ã… box and protonation) against the baseline.

    Args:
        vina (Path): Vina results directory.

    Returns:
        pd.DataFrame: One row with the paired comparison.
    """
    a = successes(load_run(vina / "baseline" / "results.csv"))
    b = successes(load_run(vina / "ligand_dimorphite_box_2" / "results.csv"))
    return pd.DataFrame([{"protocol": "box_padding_2+ligand_dimorphite", **_pair(a, b)}])


def _volume(box) -> float:
    """
    Volume of a search box.

    Args:
        box (Box): Search box.

    Returns:
        float: Volume in cubic Ã….
    """
    return float(np.prod(box.size))


def _inside(mol: Chem.Mol, box) -> bool:
    """
    Whether every heavy atom of a molecule lies inside a box.

    Args:
        mol (Chem.Mol): Molecule with a 3D conformer.
        box (Box): Search box.

    Returns:
        bool: True if no heavy atom is outside.
    """
    xyz = Chem.RemoveAllHs(mol).GetConformer().GetPositions()
    lo = np.array(box.center) - np.array(box.size) / 2
    hi = np.array(box.center) + np.array(box.size) / 2
    return bool(((xyz >= lo) & (xyz <= hi)).all())


def gyration_tables(vina: Path, dataset: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    The ligand-sized box: geometry per complex, and success against the baseline and 2 Ã… box.

    The box is rebuilt exactly as the run built it, from the crystal ligand (centre) and the
    supplied start conformer (radius of gyration).

    Args:
        vina (Path): Vina results directory.
        dataset (Path): Extracted PoseBusters root.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: Per-complex geometry (volume ratios, crystal ligand
            inside), and the comparisons: against the baseline, against the 2 Ã… box, and against
            the baseline within the complexes whose crystal ligand fits and within the rest.
    """
    rows = []
    for c in load_complexes(dataset):
        crystal = load_sdf_file(c.crystal_ligand)[0]
        start = load_sdf_file(c.start_ligand)[0]
        gyr = box_from_gyration(crystal, start)
        rows.append({
            "id": c.id,
            "ratio_to_baseline": _volume(gyr) / _volume(box_from_ligand(crystal, DEFAULT_PADDING)),
            "ratio_to_box2": _volume(gyr) / _volume(box_from_ligand(crystal, 2.0)),
            "crystal_inside": _inside(crystal, gyr),
        })  # fmt: skip
    geometry = pd.DataFrame(rows).set_index("id")
    hits = {
        r: successes(load_run(vina / r / "results.csv"))
        for r in ("baseline", "box_padding_2", "box_gyration")
    }
    gyr = hits["box_gyration"]
    comps = [
        {"comparison": "vs_baseline", "n": len(geometry), **_pair(hits["baseline"], gyr)},
        {
            "comparison": "vs_box_padding_2",
            "n": len(geometry),
            **_pair(hits["box_padding_2"], gyr),
        },
    ]
    for fits, name in ((True, "inside"), (False, "outside")):
        ids = geometry.index[geometry["crystal_inside"] == fits]
        comps.append({
            "comparison": f"vs_baseline_{name}",
            "n": len(ids),
            **_pair(hits["baseline"].loc[ids], gyr.loc[ids]),
        })  # fmt: skip
    summary = pd.DataFrame(comps)
    geometry_summary = pd.DataFrame([{
        "n": len(geometry),
        "smaller_than_baseline": int((geometry["ratio_to_baseline"] < 1).sum()),
        "median_ratio_to_baseline": geometry["ratio_to_baseline"].median(),
        "smaller_than_box2": int((geometry["ratio_to_box2"] < 1).sum()),
        "median_ratio_to_box2": geometry["ratio_to_box2"].median(),
        "crystal_partly_outside": int((~geometry["crystal_inside"]).sum()),
    }])  # fmt: skip
    return geometry_summary, summary


def main(argv: Sequence[str] | None = None) -> None:
    """
    Write the follow-up tables as CSV.

    Args:
        argv (Sequence[str] | None): Command-line arguments; sys.argv if None.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("results", type=Path, help="root with <engine>/<arm>/ runs")
    parser.add_argument("dataset", type=Path, help="extracted PoseBusters root")
    parser.add_argument("out", type=Path, help="directory for the CSV tables")
    args = parser.parse_args(argv)
    vina = args.results / "vina"
    geometry, gyration = gyration_tables(vina, args.dataset)
    tables = {
        "followup_interaction": interaction_table(vina),
        "followup_combined": combined_table(vina),
        "followup_gyration_geometry": geometry,
        "followup_gyration": gyration,
    }
    args.out.mkdir(parents=True, exist_ok=True)
    for name, table in tables.items():
        table.to_csv(args.out / f"{name}.csv", index=False)
        print(f"\n== {name}\n{table.round(4).to_string(index=False)}")


if __name__ == "__main__":
    main()
