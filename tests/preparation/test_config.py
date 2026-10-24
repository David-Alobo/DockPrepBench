import pytest

from dockprepbench.preparation.config import PrepConfig, load_config


def test_defaults_are_the_baseline():
    config = PrepConfig()

    assert (config.name, config.ph, config.box_padding) == ("baseline", 7.4, 4.0)


def test_load_config_overrides_only_given_fields(tmp_path):
    path = tmp_path / "box_padding_2.toml"
    path.write_text("box_padding = 2.0\n", encoding="utf-8")

    config = load_config(path)

    assert config == PrepConfig(name="box_padding_2", box_padding=2.0)


def test_load_config_explicit_name_wins(tmp_path):
    path = tmp_path / "file.toml"
    path.write_text('name = "acidic"\nph = 6.5\n', encoding="utf-8")

    assert load_config(path).name == "acidic"


def test_load_config_rejects_unknown_keys(tmp_path):
    path = tmp_path / "typo.toml"
    path.write_text("box_paddding = 2.0\n", encoding="utf-8")

    with pytest.raises(ValueError, match="box_paddding"):
        load_config(path)


def test_receptor_key_ignores_box_settings():
    assert PrepConfig(box_padding=2.0).receptor_key() == PrepConfig().receptor_key() == "ph7.4"
    assert PrepConfig(ph=6.5).receptor_key() == "ph6.5"


def test_shipped_configs_load():
    from pathlib import Path

    configs = sorted(Path("configs").glob("*.toml"))
    assert configs
    for path in configs:
        load_config(path)


def test_ligand_settings_default_to_supplied():
    config = PrepConfig()

    assert (config.ligand_protonation, config.ligand_conformer) == ("supplied", "supplied")


def test_invalid_ligand_setting_is_rejected(tmp_path):
    path = tmp_path / "bad.toml"
    path.write_text('ligand_conformer = "omega"\n', encoding="utf-8")

    with pytest.raises(ValueError, match="ligand_conformer"):
        load_config(path)


def test_box_size_defaults_to_padding():
    assert PrepConfig().box_size == "padding"


def test_invalid_box_size_is_rejected(tmp_path):
    path = tmp_path / "bad.toml"
    path.write_text('box_size = "sphere"\n', encoding="utf-8")

    with pytest.raises(ValueError, match="box_size"):
        load_config(path)


def test_receptor_key_ignores_box_size():
    assert PrepConfig(box_size="gyration").receptor_key() == PrepConfig().receptor_key()
