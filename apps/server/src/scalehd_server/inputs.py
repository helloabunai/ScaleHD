"""Sequencing runs in the server's data root. Browse subdirs, pair FASTQ R1/R2 files into samples.
All paths relative to the specified (``SCALEHD_DATA_ROOT``).
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path

from scalehd.simulate import true_genotype
from scalehd.structure import AlleleStructure

from .schemas import InputFolder, InputSample, InputSubfolder


class InputError(Exception):
    """A folder or file that isn't in the data folder, or doesn't exist."""


# <sample>_S<n>_L<lane>_R<read>_001, as a MiSeq names them
# lab feedback will probably change this logic as its from memory yeah !!
_MISEQ = re.compile(
    r"^(?P<sample>.+?)_S(?P<number>\d+)_L\d{3}_R(?P<read>[12])_001\.(?:fastq|fq)(?:\.gz)?$"
)
# <sample>_R1, <sample>.R1, and <sample>_1 as the names pairs.
_PLAIN = re.compile(r"^(?P<sample>.+?)[._]R?(?P<read>[12])\.(?:fastq|fq)(?:\.gz)?$")
# Any other FASTQ file = a sample of its own, single-end.
_FASTQ = re.compile(r"^(?P<sample>.+?)\.(?:fastq|fq)(?:\.gz)?$")


@dataclass(frozen=True, slots=True)
class SampleFiles:
    """One sample's FASTQ files. ``skipped`` says why they can't be run."""

    name: str
    files: list[str]  # R1 first
    r1: str | None = None
    r2: str | None = None
    skipped: str | None = None
    number: int | None = None  # the MiSeq S number


def pair_files(names: list[str]) -> list[SampleFiles]:
    """FASTQ file names paired into samples.

    An R1 without an R2 is a single-end sample. An R2 without its R1, or a sample and
    read with more than one file, is skipped with stated reason.
    Samples sorted by the run's S order (if present), then alphabetically.
    """
    groups: dict[tuple[str, int | None], dict[int, list[str]]] = {}
    for name in sorted(names):
        if match := _MISEQ.match(name):
            key: tuple[str, int | None] = (match["sample"], int(match["number"]))
            read = int(match["read"])
        elif match := _PLAIN.match(name):
            key, read = (match["sample"], None), int(match["read"])
        elif match := _FASTQ.match(name):
            key, read = (match["sample"], None), 1
        else:
            continue
        groups.setdefault(key, {}).setdefault(read, []).append(name)

    samples: list[SampleFiles] = []
    for (sample, number), reads in groups.items():
        r1s, r2s = reads.get(1, []), reads.get(2, [])
        files = r1s + r2s
        if len(r1s) > 1 or len(r2s) > 1:
            reason = "more than one file for the same sample name"
            samples.append(SampleFiles(sample, files, skipped=reason, number=number))
        elif not r1s:
            reason = "R2 without its R1"
            samples.append(SampleFiles(sample, files, skipped=reason, number=number))
        else:
            r2 = r2s[0] if r2s else None
            samples.append(SampleFiles(sample, files, r1s[0], r2, number=number))
    samples.sort(key=lambda p: (p.number is None, p.number or 0, p.name))
    return samples


def resolve_folder(root: Path, folder: str) -> Path:
    """``folder`` (relative to the data folder) as a real path, if it's a folder in it."""
    root = root.resolve()
    path = (root / folder).resolve()
    if not path.is_relative_to(root) or not path.is_dir():
        raise InputError(f"no folder {folder!r} in the data folder")
    return path


def resolve_file(root: Path, file: str) -> Path:
    """``file`` (relative to the data folder) as a real path, if it's a file in it."""
    root = root.resolve()
    path = (root / file).resolve()
    if not path.is_relative_to(root) or not path.is_file():
        raise InputError(f"no file {file!r} in the data folder")
    return path


def read_truth(r1: Path, sample: str) -> str | None:
    """The genotype a simulated sample was made from, if ``<sample>.truth.json`` is beside it.

    ``scalehd simulate`` and ``simulate-run`` write one. just for development work, real
    data won't have these.
    """
    path = r1.parent / f"{sample}.truth.json"
    try:
        data = json.loads(path.read_text())
        if data.get("schema") != "scalehd.simulation/1":
            return None
        structures = [AlleleStructure.from_label(a["structure"]) for a in data["alleles"]]
    except (OSError, ValueError, KeyError, TypeError):
        return None
    return "/".join(structure.label for structure in true_genotype(structures))


def _entries(root: Path, path: Path) -> tuple[list[str], list[str]]:
    """The folders and files in ``path``, by name. ignores hidden + symlinks."""
    folders: list[str] = []
    files: list[str] = []
    for entry in sorted(path.iterdir()):
        if entry.name.startswith("."):
            continue
        target = entry.resolve()
        if not target.is_relative_to(root):
            continue
        if target.is_dir():
            folders.append(entry.name)
        elif target.is_file():
            files.append(entry.name)
    return folders, files


def _subfolder(root: Path, path: Path) -> InputSubfolder:
    try:
        folders, files = _entries(root, path)
    except OSError:  # can't read the dir for w/e reason
        folders, files = [], []
    samples = [s for s in pair_files(files) if not s.skipped and s.name != "Undetermined"]
    return InputSubfolder(name=path.name, folders=len(folders), samples=len(samples))


def list_folder(root: Path, folder: str) -> InputFolder:
    """One folder of the data folder/it's subfolders and FASTQ samples."""
    root = root.resolve()
    path = resolve_folder(root, folder)
    folders, files = _entries(root, path)

    found = pair_files(files)
    used = {name for sample in found for name in sample.files}
    relative = path.relative_to(root).as_posix()
    prefix = "" if relative == "." else f"{relative}/"
    samples = [
        InputSample(
            name=sample.name,
            files=sample.files,
            r1=prefix + sample.r1 if sample.r1 else None,
            r2=prefix + sample.r2 if sample.r2 else None,
            size=sum((path / f).stat().st_size for f in sample.files),
            undetermined=sample.name == "Undetermined",
            skipped=sample.skipped,
        )
        for sample in found
    ]
    return InputFolder(
        folder="" if relative == "." else relative,
        path=str(path),
        folders=[_subfolder(root, path / name) for name in folders],
        samples=samples,
        other_files=[f for f in files if f not in used],
    )
