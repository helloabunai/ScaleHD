"""Users!! Register, log in and out, and change password."""

from fastapi import APIRouter, HTTPException, Request, Response, status
from sqlalchemy import func, select
from sqlalchemy.exc import IntegrityError

from ..auth import (
    CurrentUser,
    check_password,
    end_other_sessions,
    end_session,
    hash_password,
    start_session,
)
from ..config import ServerConfig
from ..db import DbSession
from ..models import User
from ..schemas import Login, NewAccount, PasswordChange, Registration, UserOut

router = APIRouter(prefix="/auth", tags=["accounts"])


@router.get("/registration")
def registration(session: DbSession, config: ServerConfig) -> Registration:
    first = session.scalar(select(func.count()).select_from(User)) == 0
    return Registration(open=first or config.allow_registration, first_account=first)


@router.post("/register", status_code=status.HTTP_201_CREATED)
def register(
    account: NewAccount,
    request: Request,
    response: Response,
    session: DbSession,
    config: ServerConfig,
) -> UserOut:
    """Create an account and log in with it. The first account = scalehd server admin."""
    first = session.scalar(select(func.count()).select_from(User)) == 0
    if not (first or config.allow_registration):
        raise HTTPException(status.HTTP_403_FORBIDDEN, "registration is closed")
    user = User(
        username=account.username, password_hash=hash_password(account.password), is_admin=first
    )
    session.add(user)
    try:
        session.flush()
    except IntegrityError:
        raise HTTPException(status.HTTP_409_CONFLICT, "that username is taken") from None
    start_session(session, user, request, response, config.session_days)
    session.commit()
    return UserOut.model_validate(user)


@router.post("/login")
def login(
    credentials: Login,
    request: Request,
    response: Response,
    session: DbSession,
    config: ServerConfig,
) -> UserOut:
    user = check_password(session, credentials.username, credentials.password)
    if user is None:
        raise HTTPException(status.HTTP_401_UNAUTHORIZED, "wrong username or password")
    start_session(session, user, request, response, config.session_days)
    session.commit()
    return UserOut.model_validate(user)


@router.post("/logout", status_code=status.HTTP_204_NO_CONTENT)
def logout(request: Request, response: Response, session: DbSession) -> None:
    end_session(session, request, response)
    session.commit()


@router.get("/me")
def me(user: CurrentUser) -> UserOut:
    return UserOut.model_validate(user)


@router.put("/password", status_code=status.HTTP_204_NO_CONTENT)
def change_password(
    change: PasswordChange, user: CurrentUser, request: Request, session: DbSession
) -> None:
    """Change the password, and log out every other browser presenting this token/user."""
    if check_password(session, user.username, change.current_password) is None:
        raise HTTPException(status.HTTP_400_BAD_REQUEST, "current password is wrong")
    user.password_hash = hash_password(change.new_password)
    end_other_sessions(session, user, request)
    session.commit()
