"""Each user's default job settings, and how new jobs pick them up."""

from fastapi.testclient import TestClient

JOB = {"name": "run 1", "samples": [{"name": "s1", "r1": "s1_R1.fastq.gz"}]}


def test_settings_need_a_login(client: TestClient) -> None:
    assert client.get("/api/settings").status_code == 401


def test_new_users_default_to_legacy_genotyping(logged_in: TestClient) -> None:
    assert logged_in.get("/api/settings").json()["method"] == "legacy"


def test_saved_settings_come_back(logged_in: TestClient) -> None:
    settings = logged_in.get("/api/settings").json() | {"method": "model", "min_molecules": 200}
    response = logged_in.put("/api/settings", json=settings)
    assert response.status_code == 200
    assert response.json() == settings
    assert logged_in.get("/api/settings").json() == settings


def test_each_user_has_their_own_settings(logged_in: TestClient) -> None:
    logged_in.put("/api/settings", json={"method": "model"})
    bob = TestClient(logged_in.app)
    bob.post("/api/auth/register", json={"username": "bob", "password": "correct horse"})
    assert bob.get("/api/settings").json()["method"] == "legacy"
    assert logged_in.get("/api/settings").json()["method"] == "model"


def test_unknown_method_is_rejected(logged_in: TestClient) -> None:
    assert logged_in.put("/api/settings", json={"method": "magic"}).status_code == 422


def test_jobs_use_the_users_default_method(logged_in: TestClient) -> None:
    response = logged_in.post("/api/jobs", json=JOB)
    assert response.status_code == 400
    assert response.json()["detail"].startswith("legacy genotyping is not available yet")

    logged_in.put("/api/settings", json={"method": "model"})
    # Past the method check; the rest of job creation isn't built yet.
    assert logged_in.post("/api/jobs", json=JOB).status_code == 501


def test_a_job_can_choose_its_own_method(logged_in: TestClient) -> None:
    assert (
        logged_in.post("/api/jobs", json=JOB | {"settings": {"method": "model"}}).status_code == 501
    )
    logged_in.put("/api/settings", json={"method": "model"})
    response = logged_in.post("/api/jobs", json=JOB | {"settings": {"method": "legacy"}})
    assert response.status_code == 400
