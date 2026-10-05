"""HTTP routes, all mounted under ``/api``."""

from fastapi import APIRouter

from . import accounts, admin, folders, health, inputs, jobs, settings, tags

api = APIRouter(prefix="/api")
for _module in (health, accounts, admin, folders, inputs, jobs, settings, tags):
    api.include_router(_module.router)
