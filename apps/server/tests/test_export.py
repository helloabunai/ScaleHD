"""A job's results as a CSV table."""

import csv
import io

from scalehd.counts import SampleCounts
from scalehd.genotype import call_genotype
from scalehd.structure import AlleleStructure, Observation
from scalehd_server.export import HEADER, results_csv, results_filename
from scalehd_server.models import Job, Sample, SampleStatus, User
from scalehd_server.schemas import JobSettings


def _call() -> dict[str, object]:
    counts = SampleCounts()
    for _ in range(600):
        counts.add(Observation.exact(AlleleStructure(17)))
        counts.add(Observation.exact(AlleleStructure(43)))
    return call_genotype(counts).to_dict()


def _job(*samples: Sample) -> Job:
    job = Job(
        id=7,
        owner=User(username="autotest-user", password_hash="x"),
        name="Run 1: cohort A",
        settings=JobSettings().model_dump(mode="json"),
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
