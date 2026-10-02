"""The workspace i.e. output main dir. results in a folder per user, a subfolder per job.

<workspace>/<username>/<job id>-<job name>/
    job.json              id, name, owner, created, settings, samples, versions
    <sample>/counts.json
    <sample>/call.json
    <sample>/input/       simulated samples only: the generated FASTQ files
    other output data will appear here as work progresses
"""

from __future__ import annotations

import json
import re
import shutil
import unicodedata
from collections.abc import Sequence
from pathlib import Path

import scalehd

from . import __version__
from .models import Job

_NOT_LETTER_OR_DIGIT = re.compile(r"[^a-z0-9]+")
_LONGEST_NAME = 60


class WorkspaceError(Exception):
    """The workspace can't be written to for whatever reason."""


def folder_name(name: str, fallback: str) -> str:
    """generate dir for jobs"""
    ascii_name = unicodedata.normalize("NFKD", name).encode("ascii", "ignore").decode()
    slug = _NOT_LETTER_OR_DIGIT.sub("-", ascii_name.lower()).strip("-")
    return slug[:_LONGEST_NAME].strip("-") or fallback


def sample_folder_names(names: Sequence[str]) -> list[str]:
    """Folder names for a job's samples, in order, made unique with -2, -3, etc"""
    used: set[str] = set()
    folders = []
    for name in names:
        base = candidate = folder_name(name, "sample")
        n = 1
        while candidate in used:
            n += 1
            candidate = f"{base}-{n}"
        used.add(candidate)
        folders.append(candidate)
    return folders


def user_folder(workspace: Path, username: str) -> Path:
    """Where a user's (user on the docker server) jobs go."""
    return workspace / username


def write_job_folder(workspace: Path, job: Job) -> Path:
    """Create a stored job's folder and its job.json, and return the folder."""
    label = "demo" if job.demo else folder_name(job.name, "job")
    folder = user_folder(workspace, job.owner.username) / f"{job.id}-{label}"
    names = sample_folder_names([sample.name for sample in job.samples])
    record = {
        "schema": "scalehd.job/1",
        "id": job.id,
        "name": job.name,
        "owner": job.owner.username,
        "demo": job.demo,
        "created_at": job.created_at.isoformat(),
        "settings": job.settings,
        "samples": [
            {
                "name": sample.name,
                "folder": name,
                "r1": sample.r1,
                "r2": sample.r2,
                "simulation": sample.simulation,
                "truth": sample.truth,
            }
            for sample, name in zip(job.samples, names, strict=True)
        ],
        "versions": {"scalehd": scalehd.__version__, "scalehd-server": __version__},
    }
    try:
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "job.json").write_text(json.dumps(record, indent=2) + "\n")
    except OSError as exc:
        reason = exc.strerror or str(exc)
        raise WorkspaceError(f"can't write to the workspace at {folder}: {reason}") from exc
    return folder


def sample_folders(job: Job) -> dict[int, Path]:
    """Each sample's folder inside the job's folder, by sample id."""
    if job.output_dir is None:
        raise ValueError(f"job {job.id} has no folder")
    names = sample_folder_names([sample.name for sample in job.samples])
    return {
        sample.id: Path(job.output_dir) / name
        for sample, name in zip(job.samples, names, strict=True)
    }


def remove_job_folder(workspace: Path, username: str, folder: Path) -> None:
    """Delete a job's folder, but only if it is inside the user's own folder.
    """
    mine = user_folder(workspace, username).resolve()
    target = folder.resolve()
    if target == mine or not target.is_relative_to(mine):
        raise WorkspaceError(f"won't delete {folder}: it is not inside {mine}")
    try:
        shutil.rmtree(target)
    except FileNotFoundError:
        return
    except OSError as exc:
        raise WorkspaceError(f"can't delete {folder}: {exc.strerror or exc}") from exc
