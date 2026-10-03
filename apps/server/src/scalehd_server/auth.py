"""Passwords and login sessions.

Passwords are stored as argon2 hashes. Logging in creates a random session token,
which goes to the browser in an HttpOnly, SameSite=Lax cookie. The database keeps only
the token's SHA-256, so logging out or expiring a session works on the server, and a
copy of the database holds no usable logins.
"""

from __future__ import annotations

import hashlib
import secrets
from datetime import UTC, datetime, timedelta
from typing import Annotated

from fastapi import Depends, HTTPException, Request, Response, status
from pwdlib import PasswordHash
from sqlalchemy import delete, select
from sqlalchemy.orm import Session

from .db import DbSession
from .models import LoginSession, User

COOKIE = "scalehd_session"

_passwords = PasswordHash.recommended()
# Checked against when a username doesn't exist, so a wrong username takes as long
# to reject as a wrong password.
_UNKNOWN_USER = _passwords.hash(secrets.token_urlsafe())


def hash_password(password: str) -> str:
    return _passwords.hash(password)


def check_password(session: Session, username: str, password: str) -> User | None:
    """The user, if the password is theirs. Rehashes it if the hash settings changed."""
    user = session.scalar(select(User).where(User.username == username))
    if user is None:
        _passwords.verify(password, _UNKNOWN_USER)
        return None
    valid, upgraded = _passwords.verify_and_update(password, user.password_hash)
    if not valid:
        return None
    if upgraded is not None:
        user.password_hash = upgraded
    return user


def start_session(
    session: Session, user: User, request: Request, response: Response, days: float
) -> None:
    """Log the user in on this browser: store a new session and set its cookie."""
    token = secrets.token_urlsafe(32)
    now = datetime.now(UTC)
    # Housekeeping i.e. this user's expired sessions.
    session.execute(
        delete(LoginSession).where(LoginSession.user_id == user.id, LoginSession.expires_at <= now)
    )
    session.add(
        LoginSession(user=user, token_hash=_digest(token), expires_at=now + timedelta(days=days))
    )
    response.set_cookie(
        COOKIE,
        token,
        max_age=int(days * 86400),
        httponly=True,
        samesite="lax",
        # idea is uni lab internal server running under a desk somewhere so Secure only
        # when the request was HTTPS.
        secure=request.url.scheme == "https",
    )


def end_session(session: Session, request: Request, response: Response) -> None:
    if token := request.cookies.get(COOKIE):
        session.execute(delete(LoginSession).where(LoginSession.token_hash == _digest(token)))
    response.delete_cookie(COOKIE, httponly=True, samesite="lax")


def end_other_sessions(session: Session, user: User, request: Request) -> None:
    """Log the user out everywhere except this browser, e.g. after a password change."""
    current = _digest(request.cookies.get(COOKIE, ""))
    session.execute(
        delete(LoginSession).where(
            LoginSession.user_id == user.id, LoginSession.token_hash != current
        )
    )


def current_user(request: Request, session: DbSession) -> User:
    """The logged-in user, from the session cookie; 401 when there is none."""
    if token := request.cookies.get(COOKIE):
        user = session.scalar(
            select(User)
            .join(LoginSession)
            .where(
                LoginSession.token_hash == _digest(token),
                LoginSession.expires_at > datetime.now(UTC),
            )
        )
        if user is not None:
            return user
    raise HTTPException(status.HTTP_401_UNAUTHORIZED, "not logged in")


def _digest(token: str) -> str:
    return hashlib.sha256(token.encode()).hexdigest()


CurrentUser = Annotated[User, Depends(current_user)]
