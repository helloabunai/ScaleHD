"""Admin view i.e. the server's users, and whoever has admin rights.
Server always has at least 1 admin (obviously).
"""

from fastapi import APIRouter, HTTPException, status
from sqlalchemy import func, select

from ..auth import AdminUser, end_all_sessions, hash_password
from ..db import DbSession
from ..models import User
from ..schemas import AdminChange, PasswordSet, UserOut

router = APIRouter(prefix="/admin", tags=["admin"])


@router.get("/users")
def list_users(admin: AdminUser, session: DbSession) -> list[UserOut]:
    users = session.scalars(select(User).order_by(User.created_at, User.id))
    return [UserOut.model_validate(user) for user in users]


@router.put("/users/{user_id}/admin")
def set_admin(user_id: int, change: AdminChange, admin: AdminUser, session: DbSession) -> UserOut:
    """Promote/demote user to/from admin. Will always be one admin."""
    user = session.get(User, user_id)
    if user is None:
        raise HTTPException(status.HTTP_404_NOT_FOUND, "no such user")
    if user.is_admin and not change.is_admin:
        admins = session.scalar(select(func.count()).select_from(User).where(User.is_admin))
        if admins == 1:
            raise HTTPException(status.HTTP_409_CONFLICT, "there must always be at least one admin")
    user.is_admin = change.is_admin
    session.commit()
    return UserOut.model_validate(user)


@router.put("/users/{user_id}/password", status_code=status.HTTP_204_NO_CONTENT)
def set_password(user_id: int, new: PasswordSet, admin: AdminUser, session: DbSession) -> None:
    user = session.get(User, user_id)
    if user is None:
        raise HTTPException(status.HTTP_404_NOT_FOUND, "no such user")
    if user.id == admin.id:
        raise HTTPException(
            status.HTTP_409_CONFLICT, "change your own password on your account page"
        )
    user.password_hash = hash_password(new.password)
    end_all_sessions(session, user)
    session.commit()
