import hashlib
import zipfile

import pytest

from dockprepbench.datasets import posebusters
from dockprepbench.datasets.posebusters import SET_NAME, download, load_complexes

SUFFIXES = ("_protein.pdb", "_ligand.sdf", "_ligand_start_conf.sdf")


def make_set(root, ids, skip=None):
    (root / SET_NAME).mkdir(parents=True)
    (root / f"{SET_NAME}_ids.txt").write_text("\n".join(ids))  # no trailing newline, as upstream
    for cid in ids:
        folder = root / SET_NAME / cid
        folder.mkdir()
        for suffix in SUFFIXES:
            if (cid, suffix) != skip:
                (folder / f"{cid}{suffix}").write_text("")


def test_load_complexes_follows_id_list(tmp_path):
    make_set(tmp_path, ["5S8I_2LY", "7PA4_C"])

    complexes = load_complexes(tmp_path)

    assert [c.id for c in complexes] == ["5S8I_2LY", "7PA4_C"]
    assert complexes[1].start_ligand.name == "7PA4_C_ligand_start_conf.sdf"


def test_load_complexes_missing_file_raises(tmp_path):
    make_set(tmp_path, ["5S8I_2LY"], skip=("5S8I_2LY", "_ligand.sdf"))

    with pytest.raises(FileNotFoundError, match="5S8I_2LY"):
        load_complexes(tmp_path)


@pytest.fixture
def fake_zenodo(tmp_path, monkeypatch):
    """Serve a tiny archive in place of the Zenodo download; returns a call counter."""
    source = tmp_path / "source"
    make_set(source, ["5S8I_2LY"])
    archive = tmp_path / "fake.zip"
    with zipfile.ZipFile(archive, "w") as zf:
        for path in source.rglob("*"):
            zf.write(path, path.relative_to(source))
    monkeypatch.setattr(posebusters, "ZENODO_MD5", hashlib.md5(archive.read_bytes()).hexdigest())

    calls = []

    def fake_retrieve(url, dest):
        calls.append(url)
        dest.write_bytes(archive.read_bytes())

    monkeypatch.setattr(posebusters.urllib.request, "urlretrieve", fake_retrieve)
    return calls


def test_download_extracts_and_skips_when_cached(tmp_path, fake_zenodo):
    data_dir = tmp_path / "data"

    root = download(data_dir)
    download(data_dir)

    assert [c.id for c in load_complexes(root)] == ["5S8I_2LY"]
    assert len(fake_zenodo) == 1


def test_download_rejects_bad_checksum(tmp_path, fake_zenodo, monkeypatch):
    monkeypatch.setattr(posebusters, "ZENODO_MD5", "0" * 32)

    with pytest.raises(ValueError, match="Checksum"):
        download(tmp_path / "data")
