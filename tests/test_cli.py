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


@pytest.fixture
def fake_pipeline(monkeypatch):
    """Replace data download and docking; returns the config names that were docked."""
    import pandas as pd

    from dockprepbench import benchmark
    from dockprepbench.datasets import posebusters

    docked = []

    def fake_run_benchmark(complexes, engine, config, work_dir, workers=1):
        docked.append(config.name)
        return pd.DataFrame(
            {
                "id": ["C0"],
                "top1_rmsd": [0.5],
                "top1_score": [-8.0],
                "top1_pb_valid": [True],
                "top1_failed_checks": [""],
                "n_poses": [9],
                "seconds": [1.0],
                "error": [""],
            }
        )

    monkeypatch.setattr(posebusters, "download", lambda data_dir: data_dir)
    monkeypatch.setattr(posebusters, "load_complexes", lambda root: ["C0"])
    monkeypatch.setattr(benchmark, "run_benchmark", fake_run_benchmark)
    return docked


def test_run_batches_configs_and_skips_finished(tmp_path, monkeypatch, fake_pipeline):
    monkeypatch.chdir(tmp_path)
    for name in ("arm_a", "arm_b"):
        (tmp_path / f"{name}.toml").write_text("box_padding = 2.0\n", encoding="utf-8")
    write_run(tmp_path / "results" / "vina" / "arm_a", [0.5], [True])  # already finished

    assert main(["run", "--config", "arm_a.toml", "arm_b.toml"]) == 0

    assert fake_pipeline == ["arm_b"]
    run_dir = tmp_path / "results" / "vina" / "arm_b"
    assert (run_dir / "results.csv").is_file()
    assert '"box_padding": 2.0' in (run_dir / "run.json").read_text(encoding="utf-8")


def test_run_force_reruns_finished(tmp_path, monkeypatch, fake_pipeline):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "arm_a.toml").write_text("", encoding="utf-8")
    write_run(tmp_path / "results" / "vina" / "arm_a", [0.5], [True])

    main(["run", "--config", "arm_a.toml", "--force"])

    assert fake_pipeline == ["arm_a"]


def test_run_rejects_out_with_several_configs(tmp_path, fake_pipeline):
    for name in ("a", "b"):
        (tmp_path / f"{name}.toml").write_text("", encoding="utf-8")

    with pytest.raises(SystemExit, match="single config"):
        main(
            [
                "run",
                "--config",
                str(tmp_path / "a.toml"),
                str(tmp_path / "b.toml"),
                "--out",
                str(tmp_path / "out"),
            ]
        )

    assert fake_pipeline == []


def test_figures_writes_table_and_plots(tmp_path, capsys):
    a = write_run(tmp_path / "baseline", [0.5, 3.0, 0.5], [True, True, True])
    b = write_run(tmp_path / "box_padding_2", [0.5, 0.5, 3.0], [True, False, True])

    assert main(["figures", str(a), str(b), "--out", str(tmp_path / "figs")]) == 0

    names = sorted(p.name for p in (tmp_path / "figs").iterdir())
    assert names == [
        "rmsd_cdf.pdf",
        "rmsd_cdf.png",
        "success_rates.pdf",
        "success_rates.png",
        "summary.csv",
    ]
    summary = (tmp_path / "figs" / "summary.csv").read_text(encoding="utf-8")
    assert "baseline" in summary and "box_padding_2" in summary


def test_compare_several_arms_reports_holm(tmp_path, capsys):
    base = write_run(tmp_path / "base", [0.5, 3.0, 0.5], [True, True, True])
    arm1 = write_run(tmp_path / "arm1", [0.5, 0.5, 3.0], [True, True, True])
    arm2 = write_run(tmp_path / "arm2", [0.5, 3.0, 0.5], [True, True, True])

    assert main(["compare", str(base), str(arm1), str(arm2)]) == 0

    out = capsys.readouterr().out
    assert out.count("Holm-adjusted p") == 2
    assert "(2 comparisons)" in out
