"""Admin tests for users, and who else is an admin"""

from collections.abc import Callable

from fastapi.testclient import TestClient


def _users(client: TestClient) -> dict[str, dict[str, object]]:
    return {u["username"]: u for u in client.get("/api/admin/users").json()}


def test_only_an_admin_lists_users(
    queued: TestClient, login_as: Callable[[TestClient, str], None]
) -> None:
    users = queued.get("/api/admin/users").json()
    assert [(u["username"], u["is_admin"]) for u in users] == [("autotest-user", True)]
    login_as(queued, "second-user")
    response = queued.get("/api/admin/users")
    assert response.status_code == 403
    assert response.json()["detail"] == "only an admin can do that"


def test_an_admin_can_promote_another_user(
    queued: TestClient, login_as: Callable[[TestClient, str], None]
) -> None:
    login_as(queued, "second-user")
    login_as(queued, "autotest-user")
    second = _users(queued)["second-user"]
    response = queued.put(f"/api/admin/users/{second['id']}/admin", json={"is_admin": True})
    assert response.status_code == 200
    assert response.json()["is_admin"] is True

    login_as(queued, "second-user")
    assert queued.get("/api/auth/me").json()["is_admin"] is True
    assert queued.get("/api/admin/users").status_code == 200


def test_the_last_admin_cannot_step_down(queued: TestClient) -> None:
    me = queued.get("/api/auth/me").json()
    response = queued.put(f"/api/admin/users/{me['id']}/admin", json={"is_admin": False})
    assert response.status_code == 409
    assert response.json()["detail"] == "there must always be at least one admin"
    assert queued.get("/api/auth/me").json()["is_admin"] is True


def test_an_admin_can_step_down_once_there_is_another(
    queued: TestClient, login_as: Callable[[TestClient, str], None]
) -> None:
    login_as(queued, "second-user")
    login_as(queued, "autotest-user")
    users = _users(queued)
    queued.put(f"/api/admin/users/{users['second-user']['id']}/admin", json={"is_admin": True})
    me = users["autotest-user"]
    response = queued.put(f"/api/admin/users/{me['id']}/admin", json={"is_admin": False})
    assert response.status_code == 200
    assert queued.get("/api/admin/users").status_code == 403


def test_promoting_an_unknown_user_is_not_found(queued: TestClient) -> None:
    assert queued.put("/api/admin/users/999/admin", json={"is_admin": True}).status_code == 404
