"""Users!! Not working yet."""

from fastapi import APIRouter, Response, status

from ..auth import CurrentUser
from ..config import ServerConfig
from ..db import DbSession
from ..errors import not_implemented
from ..schemas import Credentials, UserOut

router = APIRouter(prefix="/auth", tags=["accounts"])


@router.post("/register", status_code=status.HTTP_201_CREATED)
def register(credentials: Credentials, session: DbSession, config: ServerConfig) -> UserOut:
    """Create an account. The first one becomes the admin."""
    raise not_implemented("accounts")


@router.post("/login")
def login(credentials: Credentials, response: Response, session: DbSession) -> UserOut:
    """Check the password and set the session cookie."""
    raise not_implemented("accounts")


@router.post("/logout", status_code=status.HTTP_204_NO_CONTENT)
def logout(user: CurrentUser, session: DbSession) -> None:
    raise not_implemented("accounts")


@router.get("/me")
def me(user: CurrentUser) -> UserOut:
    return UserOut.model_validate(user)
