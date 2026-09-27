"""Uvicorn entry point exposing the ASGI app."""
from deapp.api import create_app

app = create_app()
