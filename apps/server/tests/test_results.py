"""Individual sample results tests"""

from collections import Counter
from pathlib import Path
from typing import Any

import pytest
from fastapi.testclient import TestClient
from scalehd.counts import SampleCounts
from scalehd.structure import AlleleStructure, FieldStatus, Observation
from scalehd_server.models import Job, JobStatus, Sample, SampleStatus
from scalehd_server.results import cag_ccg_cells, cag_charts, ccg_distribution, reads
from scalehd_server.schemas import GenotypeMethod, JobSettings
from scalehd_server.worker import SampleTask, run_task
from sqlalchemy import update

S = AlleleStructure.from_label
LOWER_CAG = (FieldStatus.LOWER_BOUND,) + (FieldStatus.EXACT,) * 4


def _counts() -> SampleCounts:
    counts = SampleCounts(molecules=400, dropped=3, unusable=2)
    counts.complete.update(
        {
            S("16_1_1_7_2"): 10,
            S("17_1_1_7_2"): 100,
            S("42_1_1_7_2"): 30,
            S("43_1_1_7_2"): 200,
            S("17_1_1_10_2"): 5,
        }
    )
    counts.partial[Observation((83, 1, 1, 7, 2), LOWER_CAG)] = 7
    counts.read_outcomes.update({"spanning": 340, "partial": 7, "no_anchor": 50})
    return counts


def _allele(label: str, beyond: bool = False) -> dict[str, Any]:
    s = S(label)
    return {
        "structure": label,
        "beyond_read_length": beyond,
        "cag": s.cag,
        "caacag": s.caacag,
        "ccgcca": s.ccgcca,
        "ccg": s.ccg,
        "cct": s.cct,
    }


def test_two_alleles_of_one_structure_share_a_cag_chart() -> None:
    call = {"alleles": [_allele("17_1_1_7_2"), _allele("43_1_1_7_2")]}
    (chart,) = cag_charts(_counts(), call)
    assert (chart.caacag, chart.ccgcca, chart.ccg, chart.cct) == (1, 1, 7, 2)
    assert chart.alleles == ["17_1_1_7_2", "43_1_1_7_2"]
    assert chart.called == [17, 43]
    bars = {bar.cag: (bar.molecules, bar.lower_bound) for bar in chart.bars}
    assert [bar.cag for bar in chart.bars] == list(range(16, 84))  # no gaps
    assert bars[17] == (100, 0)
    assert bars[30] == (0, 0)
    assert bars[83] == (0, 7)  # only known to be at least 83 TODO: verify bounds of sequencing


def test_unconfirmed_cag_counts_are_molecules_at_that_length() -> None:
    # Read to the end of repeat tract but too near the read's own end to confirm autonomously.
    counts = _counts()
    unconfirmed = (FieldStatus.UNCONFIRMED,) + (FieldStatus.EXACT,) * 4
    counts.partial[Observation((80, 1, 1, 7, 2), unconfirmed)] = 12
    counts.partial[Observation((43, 1, 1, 7, 2), unconfirmed)] = 4
    call = {"alleles": [_allele("17_1_1_7_2"), _allele("43_1_1_7_2")]}
    (chart,) = cag_charts(counts, call)
    bars = {bar.cag: (bar.molecules, bar.lower_bound) for bar in chart.bars}
    assert bars[80] == (12, 0)
    assert bars[43] == (204, 0)
    assert bars[83] == (0, 7)


def test_alleles_on_different_ccgs_get_a_chart_each() -> None:
    call = {"alleles": [_allele("17_1_1_10_2"), _allele("43_1_1_7_2")]}
    charts = cag_charts(_counts(), call)
    assert [chart.ccg for chart in charts] == [10, 7]
    assert [bar.cag for bar in charts[0].bars] == [17]


def test_a_homozygote_has_one_chart_with_one_allele() -> None:
    call = {"alleles": [_allele("17_1_1_7_2"), _allele("17_1_1_7_2")]}
    (chart,) = cag_charts(_counts(), call)
    assert chart.alleles == ["17_1_1_7_2"]
    assert chart.called == [17]


def test_no_call_means_no_cag_charts() -> None:
    assert cag_charts(_counts(), None) == []


def test_ccg_distribution_fills_gaps() -> None:
    assert [(bar.ccg, bar.molecules) for bar in ccg_distribution(_counts())] == [
        (7, 340),
        (8, 0),
        (9, 0),
        (10, 5),
    ]


def test_cag_ccg_cells_are_the_non_empty_ones() -> None:
    cells = {(cell.cag, cell.ccg): cell.molecules for cell in cag_ccg_cells(_counts())}
    assert cells == {(16, 7): 10, (17, 7): 100, (42, 7): 30, (43, 7): 200, (17, 10): 5}


def test_reads_summarise_what_happened_to_them() -> None:
    summary = reads(_counts())
    assert (summary.molecules, summary.complete, summary.partial) == (400, 345, 7)
    assert (summary.dropped, summary.unusable) == (3, 2)
    assert summary.read_outcomes == {"spanning": 340, "partial": 7, "no_anchor": 50}


# Through the API, with one demo sample really counted and called.


def _finish_first_sample(client: TestClient) -> tuple[dict[str, Any], dict[str, Any]]:
    job = client.post("/api/jobs/demo").json()
    sample = job["samples"][0]
    folder = Path(job["output_dir"]) / sample["name"]
    result = run_task(
        SampleTask(
            sample_id=sample["id"],
            name=sample["name"],
            folder=folder,
            settings=JobSettings(method=GenotypeMethod.MODEL),
            simulation={"alleles": ["17_1_1_7_2", "21_1_1_7_2"], "pairs": 1000, "seed": 1},
            truth=sample["truth"],
        )
    )
    with client.app.state.sessions() as session:
        session.execute(
            update(Sample)
            .where(Sample.id == sample["id"])
            .values(status=SampleStatus.FINISHED, call=result.call, genotype=result.genotype)
        )
        session.execute(update(Job).where(Job.id == job["id"]).values(status=JobStatus.FINISHED))
        session.commit()
    return job, sample


def test_sample_page_data_for_a_finished_sample(quiet_client: TestClient) -> None:
    job, sample = _finish_first_sample(quiet_client)
    response = quiet_client.get(f"/api/jobs/{job['id']}/samples/{sample['id']}")
    assert response.status_code == 200
    detail = response.json()
    assert detail["sample"]["name"] == "normal-heterozygote"
    assert detail["job_name"] == job["name"]
    assert detail["call"]["genotype"] == "17_1_1_7_2/21_1_1_7_2"
    (chart,) = detail["cag_charts"]
    assert chart["called"] == [17, 21]
    assert sum(bar["molecules"] for bar in chart["bars"]) > 0
    assert any(bar["ccg"] == 7 for bar in detail["ccg"])
    assert detail["reads"]["molecules"] > 0
    assert sorted(detail["files"]) == ["call", "counts", "r1", "r2"]
    assert detail["previous_id"] is None
    assert detail["next_id"] == job["samples"][1]["id"]


def test_a_sample_that_has_not_run_has_no_distributions(quiet_client: TestClient) -> None:
    job = quiet_client.post("/api/jobs/demo").json()
    last = job["samples"][-1]
    detail = quiet_client.get(f"/api/jobs/{job['id']}/samples/{last['id']}").json()
    assert detail["call"] is None
    assert detail["cag_charts"] == detail["ccg"] == detail["cells"] == []
    assert detail["reads"] is None
    assert detail["files"] == []
    assert detail["next_id"] is None


@pytest.mark.parametrize(
    ("kind", "starts"), [("call", b"{"), ("counts", b"{"), ("r1", b"\x1f\x8b")]
)
def test_sample_files_download(quiet_client: TestClient, kind: str, starts: bytes) -> None:
    job, sample = _finish_first_sample(quiet_client)
    response = quiet_client.get(f"/api/jobs/{job['id']}/samples/{sample['id']}/files/{kind}")
    assert response.status_code == 200
    assert response.content.startswith(starts)
    assert "attachment" in response.headers["content-disposition"]


def test_missing_and_unknown_files_are_404(quiet_client: TestClient) -> None:
    job = quiet_client.post("/api/jobs/demo").json()
    sample = job["samples"][-1]
    base = f"/api/jobs/{job['id']}/samples/{sample['id']}/files"
    assert quiet_client.get(f"{base}/call").status_code == 404
    assert quiet_client.get(f"{base}/passwd").status_code == 422


def test_samples_are_private_and_belong_to_their_job(quiet_client: TestClient) -> None:
    first = quiet_client.post("/api/jobs/demo").json()
    second = quiet_client.post("/api/jobs/demo").json()
    sample = first["samples"][0]["id"]
    assert quiet_client.get(f"/api/jobs/{second['id']}/samples/{sample}").status_code == 404
    bob = TestClient(quiet_client.app)
    bob.post("/api/auth/register", json={"username": "bob", "password": "correct horse"})
    assert bob.get(f"/api/jobs/{first['id']}/samples/{sample}").status_code == 404
    assert bob.get(f"/api/jobs/{first['id']}/samples/{sample}/files/call").status_code == 404


def test_cag_chart_counts_match_the_saved_counts(quiet_client: TestClient) -> None:
    job, sample = _finish_first_sample(quiet_client)
    saved = SampleCounts.read_json(Path(job["output_dir"]) / sample["name"] / "counts.json")
    detail = quiet_client.get(f"/api/jobs/{job['id']}/samples/{sample['id']}").json()
    on_ccg7 = Counter({s.cag: n for s, n in saved.complete.items() if s.counts[1:] == (1, 1, 7, 2)})
    (chart,) = detail["cag_charts"]
    assert {bar["cag"]: bar["molecules"] for bar in chart["bars"] if bar["molecules"]} == dict(
        on_ccg7
    )
