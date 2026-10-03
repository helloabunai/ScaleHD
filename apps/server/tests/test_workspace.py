"""Folder names and the layout of a job's folder."""

import json
import os
from pathlib import Path

import pytest
import scalehd
import scalehd_server
from scalehd_server.models import Job, Sample, User
from scalehd_server.schemas import GenotypeMethod, JobSettings
from scalehd_server.workspace import (
    WorkspaceError,
    folder_name,
    remove_job_folder,
    sample_folder_names,
    sample_folders,
    write_job_folder,
)
from sqlalchemy.orm import Session, sessionmaker


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("Run 42 / MiSeq", "run-42-miseq"),
        ("", "job"),
        ("!!!", "job"),
        ("Ünïcode Bób", "unicode-bob"),
        ("../../etc/passwd", "etc-passwd"),
        ("a" * 80, "a" * 60),
        ("x" * 59 + " tail", "x" * 59),
    ],
)
def test_folder_name_is_safe(name: str, expected: str) -> None:
    assert folder_name(name, "job") == expected


def test_sample_folder_names_are_unique() -> None:
    assert sample_folder_names(["A", "a", "a-2", "!!"]) == ["a", "a-2", "a-2-2", "sample"]


def _job(session: Session, *, demo: bool, name: str = "Demo: 9 simulated samples") -> Job:
    job = Job(
        owner=User(username="autotest-user", password_hash="x"),
        name=name,
        demo=demo,
        settings=JobSettings(method=GenotypeMethod.MODEL).model_dump(mode="json"),
    )
    job.samples = [
        Sample(
            name="expanded",
            simulation={"alleles": ["17_1_1_7_2", "43_1_1_7_2"], "pairs": 5000, "seed": 2},
            truth="17_1_1_7_2/43_1_1_7_2",
        ),
        Sample(name="Expanded"),
    ]
    session.add(job)
    session.flush()
    return job


def test_job_folder_and_record(sessions: sessionmaker[Session], tmp_path: Path) -> None:
    with sessions() as session:
        job = _job(session, demo=True)
        folder = write_job_folder(tmp_path, job)
        assert folder == tmp_path / "autotest-user" / f"{job.id}-demo"
        record = json.loads((folder / "job.json").read_text())
        job.output_dir = str(folder)
        folders = sample_folders(job)

    assert record["owner"] == "autotest-user"
    assert record["settings"]["method"] == "model"
    assert [s["folder"] for s in record["samples"]] == ["expanded", "expanded-2"]
    assert record["samples"][0]["truth"] == "17_1_1_7_2/43_1_1_7_2"
    assert record["versions"] == {
        "scalehd": scalehd.__version__,
        "scalehd-server": scalehd_server.__version__,
    }
    assert sorted(folders.values()) == [folder / "expanded", folder / "expanded-2"]


def test_a_real_job_folder_uses_the_job_name(
    sessions: sessionmaker[Session], tmp_path: Path
) -> None:
    with sessions() as session:
        job = _job(session, demo=False, name="Run 42")
        assert write_job_folder(tmp_path, job).name == f"{job.id}-run-42"


@pytest.mark.skipif(os.geteuid() == 0, reason="root can write anywhere")
def test_unwritable_workspace_names_the_folder(
    sessions: sessionmaker[Session], tmp_path: Path
) -> None:
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    workspace.chmod(0o500)
    try:
        with sessions() as session, pytest.raises(WorkspaceError) as caught:
            write_job_folder(workspace, _job(session, demo=True))
    finally:
        workspace.chmod(0o700)
    assert str(caught.value).startswith(
        f"can't write to the workspace at {workspace}/autotest-user/"
    )


def test_a_job_folder_that_links_outside_your_folder_is_not_deleted(tmp_path: Path) -> None:
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    mine = tmp_path / "workspace" / "autotest-user"
    mine.mkdir(parents=True)
    (mine / "1-demo").symlink_to(elsewhere)
    with pytest.raises(WorkspaceError, match="won't delete"):
        remove_job_folder(tmp_path / "workspace", "autotest-user", mine / "1-demo")
    assert elsewhere.exists()


def test_your_own_folder_itself_is_not_deleted(tmp_path: Path) -> None:
    mine = tmp_path / "workspace" / "autotest-user"
    mine.mkdir(parents=True)
    with pytest.raises(WorkspaceError, match="won't delete"):
        remove_job_folder(tmp_path / "workspace", "autotest-user", mine)
    assert mine.exists()
