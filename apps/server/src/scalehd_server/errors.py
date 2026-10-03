"""Errors shared by the routes."""

from fastapi import HTTPException, status


def not_implemented(feature: str) -> HTTPException:
    """For routes that are declared but not built yet."""
    return HTTPException(status.HTTP_501_NOT_IMPLEMENTED, f"{feature}: not implemented yet")
