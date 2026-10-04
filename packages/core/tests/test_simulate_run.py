"""Test the simulated "real job" run"""

import json
from pathlib import Path

import pytest
from scalehd.cli import main
from scalehd.seqio import read_fastq
from scalehd.simulate_run import RunSample, placeholder_samples, write_run
from scalehd.structure import AlleleStructure


def test_placeholder_samples_cover_a_wide_range() -> None:
    samples = placeholder_samples(random_count=20, seed=1)
    names = [s.name for s in samples]
    assert len(names) == len(set(names))
    structures = [AlleleStructure.from_label(label) for s in samples for label in s.alleles]
    cags = {s.cag for s in structures}
    assert min(cags) <= 12  # short normal
    assert max(cags) >= 100  # beyond read length
    assert any(27 <= c <= 35 for c in cags)  # intermediate
    assert any(36 <= c <= 39 for c in cags)  # reduced penetrance
    assert any(60 <= c <= 80 for c in cags)  # juvenile onset
    assert {s.ccg for s in structures} >= {6, 7, 9, 10, 12}
    assert {s.caacag for s in structures} >= {0, 1, 2}
    assert {s.ccgcca for s in structures} >= {0, 1, 2}
    assert any(s.cct != 2 for s in structures)
    assert any(len(s.alleles) == 1 for s in samples)  # homozygous
    assert any(s.single_end for s in samples)
    assert any(s.folder for s in samples)
    assert any(s.pairs is not None and s.pairs < 1000 for s in samples)  # low depth
    assert placeholder_samples(random_count=20, seed=1) == samples  # the same every time


def test_write_run_lays_a_folder_out_as_a_miseq_run(tmp_path: Path) -> None:
    samples = [
        RunSample("a-17-43", ("17_1_1_7_2", "43_1_1_7_2")),
        RunSample("b-r1-only", ("18_1_1_7_2",), single_end=True),
        RunSample("c-rerun", ("20_1_1_7_2", "40_1_1_7_2"), folder="reruns"),
    ]
    run = tmp_path / "run-01"
    write_run(run, samples, pairs=40, seed=3, workers=1)

    files = sorted(p.relative_to(run).as_posix() for p in run.rglob("*") if p.is_file())
    assert files == sorted(
        [
            "README.txt",
            "Undetermined_S0_L001_R1_001.fastq.gz",
            "Undetermined_S0_L001_R2_001.fastq.gz",
            "a-17-43_S1_L001_R1_001.fastq.gz",
            "a-17-43_S1_L001_R2_001.fastq.gz",
            "a-17-43.truth.json",
            "b-r1-only_S2_L001_R1_001.fastq.gz",
            "b-r1-only.truth.json",
            "reruns/c-rerun_S3_L001_R1_001.fastq.gz",
            "reruns/c-rerun_S3_L001_R2_001.fastq.gz",
            "reruns/c-rerun.truth.json",
            # An R2 whose R1 is missing (which we don't support)
            "orphan_S4_L001_R2_001.fastq.gz",
        ]
    )
    truth = json.loads((run / "a-17-43.truth.json").read_text())
    assert [a["structure"] for a in truth["alleles"]] == ["17_1_1_7_2", "43_1_1_7_2"]
    assert len(list(read_fastq(run / "a-17-43_S1_L001_R1_001.fastq.gz"))) == 40
    assert "generated" in (run / "README.txt").read_text().lower()


def test_write_run_refuses_a_folder_that_has_files(tmp_path: Path) -> None:
    (tmp_path / "keep.txt").write_text("")
    with pytest.raises(FileExistsError):
        write_run(tmp_path, [RunSample("a", ("17_1_1_7_2",))], pairs=10, workers=1)
    assert [p.name for p in tmp_path.iterdir()] == ["keep.txt"]


def test_cli_simulate_run(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    run = tmp_path / "run-01"
    assert main(["simulate-run", str(run), "--pairs", "20", "--random", "2", "--workers", "1"]) == 0
    assert (run / "README.txt").is_file()
    assert len(list(run.rglob("*.truth.json"))) == len(placeholder_samples(random_count=2))
    assert str(run) in capsys.readouterr().out
