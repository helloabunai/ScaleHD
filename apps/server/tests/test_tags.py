"""Test tags!! shared labels for jobs, e.g. the paper or cohort they're for."""

from collections.abc import Callable
from typing import Any

from fastapi.testclient import TestClient

SAMPLE = {"name": "a", "r1": "run-01/a_R1.fastq.gz", "r2": "run-01/a_R2.fastq.gz"}


def _tag(client: TestClient, name: str) -> dict[str, Any]:
    response = client.post("/api/tags", json={"name": name})
    assert response.status_code in (200, 201), response.text
    tag: dict[str, Any] = response.json()
    return tag


def _job(client: TestClient, tags: list[int]) -> Any:
    body = {"name": "run", "samples": [SAMPLE], "settings": {"method": "model"}, "tags": tags}
    return client.post("/api/jobs", json=body)


def test_tags_need_a_login(client: TestClient) -> None:
    assert client.get("/api/tags").status_code == 401
    assert client.post("/api/tags", json={"name": "x"}).status_code == 401


def test_anyone_can_make_a_tag_and_everyone_sees_it(
    queued: TestClient, login_as: Callable[[TestClient, str], None]
) -> None:
    login_as(queued, "second-user")
    response = queued.post("/api/tags", json={"name": "  Paper XYZ  "})
    assert response.status_code == 201
    tag = response.json()
    assert tag["name"] == "Paper XYZ"
    assert tag["jobs"] == 0

    login_as(queued, "autotest-user")
    assert queued.get("/api/tags").json() == [tag]
    # The same name in other case is the same tag.
    again = queued.post("/api/tags", json={"name": "paper xyz"})
    assert again.status_code == 200
    assert again.json()["id"] == tag["id"]


def test_tag_names_are_1_to_15_characters(queued: TestClient) -> None:
    assert queued.post("/api/tags", json={"name": "   "}).status_code == 422
    assert queued.post("/api/tags", json={"name": "x" * 16}).status_code == 422
    assert queued.post("/api/tags", json={"name": "x" * 15}).status_code == 201


def test_a_job_takes_up_to_five_tags(queued: TestClient) -> None:
    ids = [_tag(queued, name)["id"] for name in ("e", "b", "a", "d", "c", "f")]
    response = _job(queued, ids[:5])
    assert response.status_code == 201
    assert [t["name"] for t in response.json()["tags"]] == ["a", "b", "c", "d", "e"]
    assert _job(queued, ids).status_code == 422
    unknown = _job(queued, [ids[0], 999])
    assert unknown.status_code == 400
    assert unknown.json()["detail"] == "no such tag: 999"
    counts = {t["name"]: t["jobs"] for t in queued.get("/api/tags").json()}
    assert counts == {"a": 1, "b": 1, "c": 1, "d": 1, "e": 1, "f": 0}


def test_the_jobs_list_shows_tags(queued: TestClient) -> None:
    tag = _tag(queued, "cohort 2")
    _job(queued, [tag["id"]])
    assert queued.get("/api/jobs").json()[0]["tags"] == [{"id": tag["id"], "name": "cohort 2"}]


def test_only_an_admin_renames_or_deletes_tags(
    queued: TestClient, login_as: Callable[[TestClient, str], None]
) -> None:
    tag = _tag(queued, "Paper XYZ")
    other = _tag(queued, "Paper ABC")
    job = _job(queued, [tag["id"]]).json()

    login_as(queued, "second-user")
    assert queued.patch(f"/api/tags/{tag['id']}", json={"name": "x"}).status_code == 403
    assert queued.delete(f"/api/tags/{tag['id']}").status_code == 403

    login_as(queued, "autotest-user")
    renamed = queued.patch(f"/api/tags/{tag['id']}", json={"name": "Paper XYZ 2026"})
    assert renamed.status_code == 200
    assert renamed.json()["name"] == "Paper XYZ 2026"
    assert queued.get(f"/api/jobs/{job['id']}").json()["tags"][0]["name"] == "Paper XYZ 2026"
    clash = queued.patch(f"/api/tags/{tag['id']}", json={"name": "paper abc"})
    assert clash.status_code == 409

    assert queued.delete(f"/api/tags/{tag['id']}").status_code == 204
    assert queued.get(f"/api/jobs/{job['id']}").json()["tags"] == []
    assert [t["id"] for t in queued.get("/api/tags").json()] == [other["id"]]
    assert queued.delete(f"/api/tags/{tag['id']}").status_code == 404


def test_only_an_admin_changes_a_jobs_tags(
    queued: TestClient, login_as: Callable[[TestClient, str], None]
) -> None:
    ids = [_tag(queued, name)["id"] for name in ("a", "b", "c", "d", "e", "f")]
    job = _job(queued, []).json()
    changed = queued.put(f"/api/jobs/{job['id']}/tags", json={"tags": ids[:2]})
    assert changed.status_code == 200
    assert [t["name"] for t in changed.json()["tags"]] == ["a", "b"]
    assert queued.put(f"/api/jobs/{job['id']}/tags", json={"tags": ids}).status_code == 422

    login_as(queued, "second-user")
    own = _job(queued, [ids[0]]).json()
    assert queued.put(f"/api/jobs/{own['id']}/tags", json={"tags": []}).status_code == 403


def test_a_samples_results_show_its_jobs_tags(queued: TestClient) -> None:
    tag = _tag(queued, "Paper XYZ")
    job = _job(queued, [tag["id"]]).json()
    detail = queued.get(f"/api/jobs/{job['id']}/samples/{job['samples'][0]['id']}").json()
    assert detail["tags"] == [{"id": tag["id"], "name": "Paper XYZ"}]


def test_tags_sort_alphabetically_ignoring_case(queued: TestClient) -> None:
    ids = [_tag(queued, name)["id"] for name in ("Paper XYZ", "cohort 2", "Beta")]
    assert [t["name"] for t in queued.get("/api/tags").json()] == ["Beta", "cohort 2", "Paper XYZ"]
    job = _job(queued, ids).json()
    assert [t["name"] for t in job["tags"]] == ["Beta", "cohort 2", "Paper XYZ"]
