"""Browsing the data folder for FASTQ files, paired into samples."""

from collections.abc import Iterator
from pathlib import Path

import pytest
from fastapi.testclient import TestClient
from scalehd_server.app import create_app
from scalehd_server.config import ServerSettings
from scalehd_server.inputs import pair_files

RUN = [
    # Out of order on purpose. Should be returned sorted/filtered properly.
    "b-sample_S2_L001_R1_001.fastq.gz",
    "b-sample_S2_L001_R2_001.fastq.gz",
    "a-sample_S1_L001_R1_001.fastq.gz",
    "a-sample_S1_L001_R2_001.fastq.gz",
    "r1-only_S3_L001_R1_001.fastq.gz",
    "orphan_S4_L001_R2_001.fastq.gz",
    "Undetermined_S0_L001_R1_001.fastq.gz",
    "Undetermined_S0_L001_R2_001.fastq.gz",
    "a-sample.truth.json",
    "README.txt",
    ".hidden_R1.fastq.gz",
]


@pytest.fixture
def data_root(tmp_path: Path) -> Path:
    root = tmp_path / "data root"
    run = root / "run-01"
    (run / "reruns").mkdir(parents=True)
    for name in RUN:
        (run / name).write_bytes(b"@r\nACGT\n+\nIIII\n")
    (run / "reruns" / "c_R1.fastq").write_text("")
    (run / "reruns" / "c_R2.fastq").write_text("")
    (root / ".cache").mkdir()
    return root


@pytest.fixture
def browser(tmp_path: Path, data_root: Path) -> Iterator[TestClient]:
    """Logged in, on a server whose data folder is ``data_root``."""
    settings = ServerSettings(
        database_dir=tmp_path / "db",
        workspace=tmp_path / "workspace",
        data_root=data_root,
        workers=1,
    )
    with TestClient(create_app(settings)) as client:
        body = {"username": "autotest-user", "password": "correct horse"}
        assert client.post("/api/auth/register", json=body).status_code == 201
        yield client


def test_browsing_needs_a_login(client: TestClient) -> None:
    assert client.get("/api/inputs").status_code == 401


def test_a_server_without_a_data_folder_says_so(logged_in: TestClient) -> None:
    response = logged_in.get("/api/inputs")
    assert response.status_code == 409
    assert "data folder" in response.json()["detail"]


def test_the_data_folder_lists_its_folders(browser: TestClient, data_root: Path) -> None:
    listing = browser.get("/api/inputs").json()
    assert listing["folder"] == ""
    assert listing["path"] == str(data_root)
    # Hidden ones left out. Each says what it holds, for the folder tree.
    assert listing["folders"] == [{"name": "run-01", "folders": 1, "samples": 3}]
    assert listing["samples"] == []


def test_a_run_folder_pairs_its_reads_into_samples(browser: TestClient) -> None:
    listing = browser.get("/api/inputs", params={"folder": "run-01"}).json()
    assert listing["folders"] == [{"name": "reruns", "folders": 0, "samples": 1}]
    samples = {s["name"]: s for s in listing["samples"]}
    assert list(samples) == ["Undetermined", "a-sample", "b-sample", "r1-only", "orphan"]
    assert samples["a-sample"]["files"] == [
        "a-sample_S1_L001_R1_001.fastq.gz",
        "a-sample_S1_L001_R2_001.fastq.gz",
    ]
    assert samples["a-sample"]["r1"] == "run-01/a-sample_S1_L001_R1_001.fastq.gz"
    assert samples["a-sample"]["r2"] == "run-01/a-sample_S1_L001_R2_001.fastq.gz"
    assert samples["a-sample"]["size"] > 0
    assert samples["a-sample"]["skipped"] is None
    assert samples["r1-only"]["r2"] is None
    assert samples["Undetermined"]["undetermined"] is True
    assert samples["a-sample"]["undetermined"] is False
    assert listing["other_files"] == ["README.txt", "a-sample.truth.json"]


def test_files_that_cant_be_run_are_listed_with_why(browser: TestClient) -> None:
    listing = browser.get("/api/inputs", params={"folder": "run-01"}).json()
    orphan = next(s for s in listing["samples"] if s["name"] == "orphan")
    assert orphan["files"] == ["orphan_S4_L001_R2_001.fastq.gz"]
    assert orphan["skipped"] == "R2 without its R1"
    assert orphan["r1"] is None
    assert orphan["r2"] is None


def test_a_nested_folder_pairs_plain_names(browser: TestClient) -> None:
    listing = browser.get("/api/inputs", params={"folder": "run-01/reruns"}).json()
    assert [(s["name"], s["r1"], s["r2"]) for s in listing["samples"]] == [
        ("c", "run-01/reruns/c_R1.fastq", "run-01/reruns/c_R2.fastq")
    ]


@pytest.mark.parametrize("folder", ["..", "run-01/../..", "/etc", "no-such-folder"])
def test_folders_outside_the_data_folder_or_missing_are_refused(
    browser: TestClient, folder: str
) -> None:
    assert browser.get("/api/inputs", params={"folder": folder}).status_code == 404


def test_a_symlink_out_of_the_data_folder_is_refused(
    browser: TestClient, data_root: Path, tmp_path: Path
) -> None:
    (tmp_path / "elsewhere").mkdir()
    (data_root / "escape").symlink_to(tmp_path / "elsewhere")
    folders = browser.get("/api/inputs").json()["folders"]
    assert "escape" not in [f["name"] for f in folders]
    assert browser.get("/api/inputs", params={"folder": "escape"}).status_code == 404


@pytest.mark.parametrize(
    ("names", "expected"),
    [
        (["s_R1.fastq.gz", "s_R2.fastq.gz"], [("s", "s_R1.fastq.gz", "s_R2.fastq.gz")]),
        (["s.R1.fq", "s.R2.fq"], [("s", "s.R1.fq", "s.R2.fq")]),
        (["SRR1_1.fastq.gz", "SRR1_2.fastq.gz"], [("SRR1", "SRR1_1.fastq.gz", "SRR1_2.fastq.gz")]),
        (["solo.fastq.gz"], [("solo", "solo.fastq.gz", None)]),
        (["x_S9_L001_R1_001.fq.gz"], [("x", "x_S9_L001_R1_001.fq.gz", None)]),
    ],
)
def test_pairing_rules(names: list[str], expected: list[tuple[str, str, str | None]]) -> None:
    samples = pair_files(names)
    assert [(s.name, s.r1, s.r2) for s in samples] == expected
    assert [s.skipped for s in samples] == [None] * len(expected)


def test_pairing_skips_what_it_cant_pair_but_keeps_its_place() -> None:
    lanes = ["a_S1_L001_R1_001.fastq.gz", "a_S1_L002_R1_001.fastq.gz"]
    samples = pair_files(["c_S3_L001_R2_001.fastq.gz", *lanes, "b_S2_L001_R1_001.fastq.gz"])
    assert [(s.name, s.files, s.skipped) for s in samples] == [
        ("a", lanes, "more than one file for the same sample name"),
        ("b", ["b_S2_L001_R1_001.fastq.gz"], None),
        ("c", ["c_S3_L001_R2_001.fastq.gz"], "R2 without its R1"),
    ]
    assert [s.r1 for s in samples] == [None, "b_S2_L001_R1_001.fastq.gz", None]
