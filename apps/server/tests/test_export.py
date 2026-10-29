"""A job's results as a CSV table."""

import csv
import io
import json
from pathlib import Path
from typing import Any

from scalehd.counts import SampleCounts
from scalehd.genotype import call_genotype
from scalehd.structure import AlleleStructure, Observation
from scalehd_server.export import (
    HEADER,
    LAB_COLUMNS,
    VALIDATION_HEADER,
    results_csv,
    results_filename,
    validation_csv,
    validation_filename,
)
from scalehd_server.models import Job, Sample, SampleStatus, User
from scalehd_server.schemas import JobSettings
from scalehd_server.workspace import sample_folders


def _counts() -> SampleCounts:
    counts = SampleCounts()
    for cag, n in ((16, 90), (17, 600), (18, 10), (42, 200), (43, 600), (44, 60)):
        for _ in range(n):
            counts.add(Observation.exact(AlleleStructure(cag)))
    counts.read_outcomes.update({"r1.spanning": 1500, "r1.no_anchor": 12, "r2.spanning": 1490})
    counts.discordant.update({"cag": 4})
    return counts


def _call() -> dict[str, Any]:
    return call_genotype(_counts()).to_dict()


def _job(*samples: Sample, output_dir: Path | None = None) -> Job:
    job = Job(
        id=7,
        owner=User(username="autotest-user", password_hash="x"),
        name="Run 1: cohort A",
        settings=JobSettings(min_molecules=200).model_dump(mode="json"),
        output_dir=str(output_dir) if output_dir else None,
    )
    job.samples = list(samples)
    return job


def _rows(job: Job) -> list[dict[str, str]]:
    return list(csv.DictReader(io.StringIO(results_csv(job))))


def test_a_called_sample_has_both_alleles_side_by_side() -> None:
    call = _call()
    sample = Sample(
        name="s1",
        status=SampleStatus.FINISHED,
        call=call,
        truth="17_1_1_7_2/43_1_1_7_2",
        matches_truth=True,
        r1="/data/run/s1_R1.fastq.gz",
        r2="/data/run/s1_R2.fastq.gz",
    )
    (row,) = _rows(_job(sample))
    assert list(row) == HEADER
    assert (row["Sample"], row["Status"], row["Genotype"]) == (
        "s1",
        "finished",
        "17_1_1_7_2/43_1_1_7_2",
    )
    assert row["Quality"] == f"{call['quality']:.1f}"
    assert (row["Allele 1"], row["Allele 1 CAG"], row["Allele 1 CCG"]) == ("17_1_1_7_2", "17", "7")
    assert (row["Allele 2"], row["Allele 2 CAG"]) == ("43_1_1_7_2", "43")
    assert abs(float(row["Allele 1 share"]) + float(row["Allele 2 share"]) - 1) < 0.002
    assert row["Allele 2 CAG estimate"] == ""
    assert (row["Truth"], row["Matches truth"], row["Error"]) == (
        "17_1_1_7_2/43_1_1_7_2",
        "yes",
        "",
    )
    assert row["R2"] == "/data/run/s1_R2.fastq.gz"


def test_samples_not_called_are_listed_with_their_results_blank() -> None:
    rows = _rows(
        _job(
            Sample(name="waiting", status=SampleStatus.QUEUED),
            Sample(name="broken", status=SampleStatus.FAILED, error="ValueError: no molecules"),
        )
    )
    assert [(r["Sample"], r["Status"], r["Genotype"], r["Allele 1"]) for r in rows] == [
        ("waiting", "queued", "", ""),
        ("broken", "failed", "", ""),
    ]
    assert rows[1]["Error"] == "ValueError: no molecules"


def test_an_allele_beyond_read_length_shows_its_lower_bound_and_estimate() -> None:
    call = _call()
    long = dict(call["alleles"][1])  # type: ignore[index]
    long |= {"structure": "84_1_1_7_2", "cag": 84, "beyond_read_length": True}
    long |= {"cag_estimate": [93, 89, 96]}
    call["alleles"] = [call["alleles"][0], long]  # type: ignore[index]
    (row,) = _rows(_job(Sample(name="s1", status=SampleStatus.FINISHED, call=call)))
    assert (row["Allele 2 CAG"], row["Allele 2 CAG estimate"]) == ("84+", "93 (89-96)")


def test_free_text_cant_be_read_as_a_spreadsheet_formula() -> None:
    (row,) = _rows(_job(Sample(name="=HYPERLINK(1)", status=SampleStatus.QUEUED)))
    assert row["Sample"] == "'=HYPERLINK(1)"


def test_the_file_is_named_after_the_job() -> None:
    assert results_filename(_job()) == "7-run-1-cohort-a-results.csv"


def _validation_job(tmp_path: Path) -> Job:
    """A finished sample with its counts and job.json saved, and a failed one."""
    run = "/data/run-a"
    called = Sample(
        id=1,
        name="patient-0042",
        status=SampleStatus.FINISHED,
        call=_call(),
        r1=f"{run}/patient-0042_R1.fastq.gz",
        r2=f"{run}/patient-0042_R2.fastq.gz",
    )
    failed = Sample(
        id=2,
        name="patient-0043",
        status=SampleStatus.FAILED,
        error=f"FileNotFoundError: {run}/patient-0043_R1.fastq.gz",
        r1=f"{run}/patient-0043_R1.fastq.gz",
    )
    other_run = Sample(id=3, name="patient-0044", r1="/data/run-b/patient-0044_R1.fastq.gz")
    job = _job(called, failed, other_run, output_dir=tmp_path / "7-run-1")
    folder = sample_folders(job)[1]
    folder.mkdir(parents=True)
    _counts().write_json(folder / "counts.json")
    versions = {"scalehd": "2.0.0", "scalehd-server": "2.0.1"}
    (tmp_path / "7-run-1" / "job.json").write_text(json.dumps({"versions": versions}))
    return job


def test_the_validation_export_names_no_sample_file_or_job(tmp_path: Path) -> None:
    job = _validation_job(tmp_path)
    text = validation_csv(job)
    assert "patient" not in text
    assert "data/run" not in text
    assert "cohort" not in text
    assert validation_filename(job) == "scalehd-validation-job-7.csv"
    rows = list(csv.DictReader(io.StringIO(text)))
    assert list(rows[0]) == VALIDATION_HEADER
    # Only the kind of error, and the reads' folders as batches.
    assert rows[1]["Error"] == "FileNotFoundError"
    assert [row["Batch"] for row in rows] == ["B1", "B1", "B2"]
    assert all(row[column] == "" for row in rows for column in LAB_COLUMNS)


def test_both_exports_give_each_sample_the_same_validation_id(tmp_path: Path) -> None:
    job = _validation_job(tmp_path)
    results = [(r["Sample"], r["Validation ID"]) for r in _rows(job)]
    assert results == [
        ("patient-0042", "7-001"),
        ("patient-0043", "7-002"),
        ("patient-0044", "7-003"),
    ]
    validation = csv.DictReader(io.StringIO(validation_csv(job)))
    assert [r["Validation ID"] for r in validation] == ["7-001", "7-002", "7-003"]


def test_the_validation_export_has_what_comparing_and_refitting_need(tmp_path: Path) -> None:
    job = _validation_job(tmp_path)
    row = next(csv.DictReader(io.StringIO(validation_csv(job))))
    call = job.samples[0].call
    assert call is not None
    assert row["Genotype"] == "17_1_1_7_2/43_1_1_7_2"
    assert row["Alternative 1"] == call["alternatives"][0]["genotype"]
    assert float(row["Background"]) == float(f"{call['background']:.4g}")
    stutter = call["alleles"][1]["stutter"]
    assert float(row["Allele 2 stutter contraction"]) == float(f"{stutter['contraction']:.4g}")
    # Molecules of 43_1_1_7_2 from CAG 23 to 73: 200 at 42, 600 at 43, 60 at 44.
    peak = [int(n) for n in row["Allele 2 molecules N-20 to N+30"].split()]
    assert len(peak) == 51
    assert peak[19:22] == [200, 600, 60]
    assert (row["Molecules"], row["Complete"], row["Dropped"]) == ("1560", "1560", "0")
    assert (row["R1 spanning"], row["R1 no anchor"], row["R2 spanning"]) == ("1500", "12", "1490")
    assert row["Discordant"] == "cag 4"
    assert (row["ScaleHD version"], row["Server version"]) == ("2.0.0", "2.0.1")
    assert row["Method"] == "model"
    assert "min_molecules=200" in row["Settings"]
    assert "min_posterior=default" in row["Settings"]
