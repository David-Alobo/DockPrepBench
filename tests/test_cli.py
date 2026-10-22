import pytest

from dockprepbench import __version__
from dockprepbench.cli import main


def test_main_returns_zero():
    assert main([]) == 0


def test_version_flag(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--version"])

    assert exc.value.code == 0
    assert __version__ in capsys.readouterr().out


def write_run(path, rmsds, valid):
    import pandas as pd

    path.mkdir(parents=True)
    pd.DataFrame(
        {
            "id": [f"C{i}" for i in range(len(rmsds))],
            "top1_rmsd": rmsds,
            "top1_score": -8.0,
            "top1_pb_valid": valid,
            "top1_failed_checks": "",
            "n_poses": 9,
            "seconds": 1.0,
            "error": ["" if r == r else "RuntimeError: x" for r in rmsds],
        }
    ).to_csv(path / "results.csv", index=False)
    return path


def test_summarize_prints_rates(tmp_path, capsys):
    run = write_run(tmp_path / "a", [0.5, 3.0, float("nan")], [True, True, False])

    assert main(["summarize", str(run)]) == 0

    out = capsys.readouterr().out
    assert "3 complexes, 1 failed" in out
    assert "33.3%" in out
    assert "RuntimeError 1" in out


def test_compare_prints_flipped_ids(tmp_path, capsys):
    a = write_run(tmp_path / "a", [0.5, 3.0, 0.5], [True, True, True])
    b = write_run(tmp_path / "b", [0.5, 0.5, 3.0], [True, True, True])

    assert main(["compare", str(a), str(b)]) == 0

    out = capsys.readouterr().out
    assert "only A: C2" in out
    assert "only B: C1" in out
    assert "McNemar exact p = 1" in out
