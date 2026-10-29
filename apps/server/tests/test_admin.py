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


def _register(client: TestClient, username: str) -> TestClient:
    browser = TestClient(client.app)
    body = {"username": username, "password": "correct horse"}
    assert browser.post("/api/auth/register", json=body).status_code == 201
    return browser


def _login(client: TestClient, username: str, password: str) -> int:
    body = {"username": username, "password": password}
    return TestClient(client.app).post("/api/auth/login", json=body).status_code


def test_an_admin_can_set_another_users_password(queued: TestClient) -> None:
    bob = _register(queued, "bob")
    bob_id = _users(queued)["bob"]["id"]
    response = queued.put(
        f"/api/admin/users/{bob_id}/password", json={"password": "battery staple"}
    )
    assert response.status_code == 204
    assert bob.get("/api/auth/me").status_code == 401
    assert _login(queued, "bob", "correct horse") == 401
    assert _login(queued, "bob", "battery staple") == 200
    assert queued.get("/api/auth/me").status_code == 200


def test_only_an_admin_sets_other_users_passwords(queued: TestClient) -> None:
    admin_id = queued.get("/api/auth/me").json()["id"]
    bob = _register(queued, "bob")
    response = bob.put(f"/api/admin/users/{admin_id}/password", json={"password": "battery staple"})
    assert response.status_code == 403
    assert _login(queued, "autotest-user", "correct horse") == 200


def test_an_admin_changes_their_own_password_on_their_account_page(queued: TestClient) -> None:
    me = queued.get("/api/auth/me").json()
    response = queued.put(
        f"/api/admin/users/{me['id']}/password", json={"password": "battery staple"}
    )
    assert response.status_code == 409
    assert response.json()["detail"] == "change your own password on your account page"


def test_a_set_password_needs_eight_characters(queued: TestClient) -> None:
    _register(queued, "bob")
    bob_id = _users(queued)["bob"]["id"]
    response = queued.put(f"/api/admin/users/{bob_id}/password", json={"password": "short"})
    assert response.status_code == 422
    assert _login(queued, "bob", "correct horse") == 200


def test_setting_an_unknown_users_password_is_not_found(queued: TestClient) -> None:
    response = queued.put("/api/admin/users/999/password", json={"password": "battery staple"})
    assert response.status_code == 404
