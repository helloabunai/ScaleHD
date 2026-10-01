"""HTTP routes, all mounted under ``/api``."""

from fastapi import APIRouter

from . import accounts, health, inputs, jobs, settings

api = APIRouter(prefix="/api")
for _module in (health, accounts, inputs, jobs, settings):
    api.include_router(_module.router)
