"""Accounts and login sessions. Not built yet.

Plan: argon2 password hashes (``pwdlib``), and a random session token in an HttpOnly,
SameSite=Lax cookie. Tokens are stored hashed in a sessions table, so logout and
expiry work on the server. The first account registered becomes the admin, and
``SCALEHD_ALLOW_REGISTRATION=false`` closes registration after that.
"""

from __future__ import annotations

from typing import Annotated

from fastapi import Depends, Request

from .db import DbSession
from .errors import not_implemented
from .models import User


def current_user(request: Request, session: DbSession) -> User:
    """The logged-in user, from the session cookie; 401 when there is none."""
    raise not_implemented("accounts")


CurrentUser = Annotated[User, Depends(current_user)]
