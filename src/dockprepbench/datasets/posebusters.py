"""The PoseBusters benchmark set (Buttenschoen et al., 2024)."""

import hashlib
import logging
import urllib.request
import zipfile
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)

ZENODO_URL = "https://zenodo.org/api/records/8278563/files/posebusters_paper_data.zip/content"
ZENODO_MD5 = "f004ac7c4e68317b5348497d2bb6bee6"
ARCHIVE_NAME = "posebusters_paper_data.zip"
SET_NAME = "posebusters_benchmark_set"


@dataclass(frozen=True)
class Complex:
    """One protein-ligand complex, named <PDB ID>_<CCD ID>."""

    id: str
    protein: Path
    crystal_ligand: Path
    start_ligand: Path


def download(data_dir: str | Path) -> Path:
    """
    Download and extract the PoseBusters archive, verifying its checksum.

    Skips the download if a verified archive is already present.

    Args:
        data_dir (str | Path): Directory to hold the archive and extracted files.

    Returns:
        Path: Directory containing the extracted benchmark set.

    Raises:
        ValueError: If the archive checksum does not match.
    """
    data_dir = Path(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    archive = data_dir / ARCHIVE_NAME

    if not (archive.is_file() and _md5(archive) == ZENODO_MD5):
        logger.info("Downloading PoseBusters archive to %s", archive)
        urllib.request.urlretrieve(ZENODO_URL, archive)
        if _md5(archive) != ZENODO_MD5:
            raise ValueError(f"Checksum mismatch for {archive}")

    root = data_dir / "posebusters"
    if not (root / SET_NAME).is_dir():
        with zipfile.ZipFile(archive) as zf:
            zf.extractall(root)
    return root


def load_complexes(root: str | Path) -> list[Complex]:
    """
    List the complexes of the extracted PoseBusters benchmark set.

    Args:
        root (str | Path): Directory returned by `download`.

    Returns:
        list[Complex]: Complexes in the order of the official ID list.

    Raises:
        FileNotFoundError: If the ID list or any complex file is missing.
    """
    root = Path(root)
    ids = (root / f"{SET_NAME}_ids.txt").read_text().split()
    complexes = []
    for cid in ids:
        folder = root / SET_NAME / cid
        entry = Complex(
            id=cid,
            protein=folder / f"{cid}_protein.pdb",
            crystal_ligand=folder / f"{cid}_ligand.sdf",
            start_ligand=folder / f"{cid}_ligand_start_conf.sdf",
        )
        for path in (entry.protein, entry.crystal_ligand, entry.start_ligand):
            if not path.is_file():
                raise FileNotFoundError(f"Missing file for {cid}: {path}")
        complexes.append(entry)
    return complexes


def _md5(path: Path) -> str:
    """
    Compute the MD5 hash of a file.

    Args:
        path (Path): Path to the file.

    Returns:
        str: The MD5 hash of the file.
    """
    digest = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()
