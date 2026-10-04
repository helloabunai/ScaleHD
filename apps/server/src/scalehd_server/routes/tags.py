"""Tags are basic labels shared by every user, for visually grouping jobs, e.g. the paper
they're for.

Anyone can list tags and make one. only an admin can rename or delete one.
"""

from fastapi import APIRouter, HTTPException, Response, status
from sqlalchemy import func, select
from sqlalchemy.orm import selectinload

from ..auth import AdminUser, CurrentUser
from ..db import DbSession
from ..models import Tag
from ..schemas import TagChange, TagOut

router = APIRouter(prefix="/tags", tags=["tags"])


def _out(tag: Tag) -> TagOut:
    return TagOut(id=tag.id, name=tag.name, jobs=len(tag.jobs))


def _named(session: DbSession, name: str) -> Tag | None:
    """The tag with this name, ignoring case."""
    return session.scalar(select(Tag).where(func.lower(Tag.name) == name.lower()))


@router.get("")
def list_tags(user: CurrentUser, session: DbSession) -> list[TagOut]:
    tags = session.scalars(
        select(Tag).options(selectinload(Tag.jobs)).order_by(func.lower(Tag.name))
    )
    return [_out(tag) for tag in tags]


@router.post("", status_code=status.HTTP_201_CREATED)
def create_tag(
    change: TagChange, user: CurrentUser, session: DbSession, response: Response
) -> TagOut:
    """Make a tag. Return if already existing"""
    if (existing := _named(session, change.name)) is not None:
        response.status_code = status.HTTP_200_OK
        return _out(existing)
    tag = Tag(name=change.name, created_by_id=user.id)
    session.add(tag)
    session.commit()
    return _out(tag)


def _tag(session: DbSession, tag_id: int) -> Tag:
    tag = session.get(Tag, tag_id)
    if tag is None:
        raise HTTPException(status.HTTP_404_NOT_FOUND, "no such tag")
    return tag


@router.patch("/{tag_id}")
def rename_tag(tag_id: int, change: TagChange, admin: AdminUser, session: DbSession) -> TagOut:
    tag = _tag(session, tag_id)
    other = _named(session, change.name)
    if other is not None and other.id != tag.id:
        raise HTTPException(status.HTTP_409_CONFLICT, f"a tag named {other.name!r} already exists")
    tag.name = change.name
    session.commit()
    return _out(tag)


@router.delete("/{tag_id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_tag(tag_id: int, admin: AdminUser, session: DbSession) -> None:
    """Delete a tag. Jobs that had it keep running, just without the tag."""
    session.delete(_tag(session, tag_id))
    session.commit()
