"""Accounts. Registering, logging in and out, sessions, and changing password."""

from datetime import UTC, datetime, timedelta
from pathlib import Path

import pytest
from fastapi.testclient import TestClient
from httpx2 import Response
from scalehd_server.app import create_app
from scalehd_server.auth import COOKIE
from scalehd_server.config import ServerSettings
from scalehd_server.models import LoginSession, User
from sqlalchemy import select, update

PASSWORD = "correct horse"


def register(
    client: TestClient, username: str = "autotest-user", password: str = PASSWORD
) -> Response:
    return client.post("/api/auth/register", json={"username": username, "password": password})


def login(
    client: TestClient, username: str = "autotest-user", password: str = PASSWORD
) -> Response:
    return client.post("/api/auth/login", json={"username": username, "password": password})


def test_first_account_is_the_admin_and_is_logged_in(client: TestClient) -> None:
    assert client.get("/api/auth/registration").json() == {"open": True, "first_account": True}
    response = register(client)
    assert response.status_code == 201
    first = response.json()
    assert first["is_admin"]
    assert datetime.fromisoformat(first["created_at"]).tzinfo is not None
    assert client.get("/api/auth/me").json() == first
    assert client.get("/api/auth/registration").json() == {"open": True, "first_account": False}

    client.cookies.clear()
    assert not register(client, "bob").json()["is_admin"]


def test_session_cookie_is_http_only(client: TestClient) -> None:
    cookie = register(client).headers["set-cookie"].lower()
    assert "httponly" in cookie
    assert "samesite=lax" in cookie
    assert "secure" not in cookie  # plain HTTP here


def test_password_is_stored_hashed(client: TestClient) -> None:
    register(client)
    with client.app.state.sessions() as session:
        stored = session.scalars(select(User.password_hash)).one()
    assert stored.startswith("$argon2id$")
    assert PASSWORD not in stored


def test_usernames_are_case_insensitive_and_unique(client: TestClient) -> None:
    assert register(client, "Autotest-User").json()["username"] == "autotest-user"
    client.cookies.clear()
    assert register(client, "autotest-user").status_code == 409
    response = login(client, " AUTOTEST-USER ")
    assert response.status_code == 200
    assert response.json()["username"] == "autotest-user"


@pytest.mark.parametrize(
    ("username", "password"),
    [
        ("", PASSWORD),
        ("a b", PASSWORD),
        ("-autotest-user", PASSWORD),
        ("a" * 65, PASSWORD),
        ("bob", "short"),
    ],
)
def test_invalid_new_accounts_are_rejected(
    client: TestClient, username: str, password: str
) -> None:
    assert register(client, username, password).status_code == 422


def test_wrong_password_and_unknown_user_look_the_same(client: TestClient) -> None:
    register(client)
    client.cookies.clear()
    wrong_password = login(client, password="not the password")
    unknown_user = login(client, "nobody")
    assert wrong_password.status_code == unknown_user.status_code == 401
    assert wrong_password.json() == unknown_user.json()
    assert client.get("/api/auth/me").status_code == 401


def test_logout_ends_the_session_on_the_server(client: TestClient) -> None:
    register(client)
    token = client.cookies[COOKIE]
    assert client.post("/api/auth/logout").status_code == 204
    assert client.get("/api/auth/me").status_code == 401
    client.cookies.set(COOKIE, token)  # a copy of the old cookie no longer works
    assert client.get("/api/auth/me").status_code == 401


def test_expired_session_is_rejected(client: TestClient) -> None:
    register(client)
    with client.app.state.sessions() as session:
        an_hour_ago = datetime.now(UTC) - timedelta(hours=1)
        session.execute(update(LoginSession).values(expires_at=an_hour_ago))
        session.commit()
    assert client.get("/api/auth/me").status_code == 401


def test_registration_can_be_closed_after_the_admin(tmp_path: Path) -> None:
    settings = ServerSettings(
        database_dir=tmp_path / "data",
        workspace=tmp_path / "ws",
        workers=1,
        allow_registration=False,
    )
    with TestClient(create_app(settings)) as client:
        assert client.get("/api/auth/registration").json()["open"]
        assert register(client).status_code == 201
        client.cookies.clear()
        assert client.get("/api/auth/registration").json() == {
            "open": False,
            "first_account": False,
        }
        assert register(client, "bob").status_code == 403


def test_password_change_logs_out_other_browsers(client: TestClient) -> None:
    register(client)
    other_browser = TestClient(client.app)
    assert login(other_browser).status_code == 200

    change = {"current_password": "not it", "new_password": "battery staple"}
    assert client.put("/api/auth/password", json=change).status_code == 400
    change["current_password"] = PASSWORD
    assert client.put("/api/auth/password", json=change).status_code == 204

    assert client.get("/api/auth/me").status_code == 200
    assert other_browser.get("/api/auth/me").status_code == 401
    assert login(other_browser).status_code == 401
    assert login(other_browser, password="battery staple").status_code == 200
