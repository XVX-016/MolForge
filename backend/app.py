"""
MolForge Backend API
FastAPI application entrypoint
"""
import logging
from importlib import import_module
from typing import Optional

from backend.config import settings

from fastapi import FastAPI, HTTPException, Request
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from backend.db import init_db

logger = logging.getLogger(__name__)

app = FastAPI(
    title=settings.API_TITLE,
    version=settings.API_VERSION,
    description=settings.API_DESCRIPTION
)

# CORS middleware must be added before routes
# Configure CORS - use regex if available, otherwise use explicit origins
cors_kwargs = {
    "allow_credentials": True,
    "allow_methods": ["GET", "POST", "PUT", "DELETE", "OPTIONS", "PATCH"],
    "allow_headers": ["*"],
    "expose_headers": ["*"],
    "max_age": 3600,
}

if settings.CORS_ALLOW_ORIGIN_REGEX:
    cors_kwargs["allow_origin_regex"] = settings.CORS_ALLOW_ORIGIN_REGEX
    # Still include explicit origins as fallback
    if settings.CORS_ORIGINS:
        cors_kwargs["allow_origins"] = settings.CORS_ORIGINS
else:
    cors_kwargs["allow_origins"] = settings.CORS_ORIGINS if settings.CORS_ORIGINS else ["http://localhost:5173", "http://127.0.0.1:5173"]

app.add_middleware(CORSMiddleware, **cors_kwargs)

# Global exception handler to ensure CORS headers on errors
@app.exception_handler(Exception)
async def global_exception_handler(request: Request, exc: Exception):
    """Handle all exceptions and ensure CORS headers are included"""
    import traceback
    error_detail = str(exc)
    if settings.LOG_LEVEL == "DEBUG":
        error_detail += f"\n{traceback.format_exc()}"
    
    # Get origin from request headers to use in CORS response
    origin = request.headers.get("origin", "*")
    
    return JSONResponse(
        status_code=500,
        content={
            "detail": error_detail,
            "path": str(request.url),
        },
        headers={
            "Access-Control-Allow-Origin": origin,
            "Access-Control-Allow-Credentials": "true",
        }
    )


@app.exception_handler(RequestValidationError)
async def validation_exception_handler(request: Request, exc: RequestValidationError):
    """Ensure CORS headers on validation errors (422)"""
    origin = request.headers.get("origin", "*")
    return JSONResponse(
        status_code=422,
        content={"detail": exc.errors(), "body": exc.body},
        headers={
            "Access-Control-Allow-Origin": origin,
            "Access-Control-Allow-Credentials": "true",
        }
    )

@app.exception_handler(HTTPException)
async def http_exception_handler(request: Request, exc: HTTPException):
    """Ensure CORS headers on HTTP exceptions (400, 404, etc.)"""
    origin = request.headers.get("origin", "*")
    return JSONResponse(
        status_code=exc.status_code,
        content={"detail": exc.detail},
        headers={
            "Access-Control-Allow-Origin": origin,
            "Access-Control-Allow-Credentials": "true",
        }
    )

# Initialize database
init_db()

ROUTER_SPECS = [
    ("backend.routes.generate", "/generate", ["generate"]),
    ("backend.routes.library", None, ["molecules"]),
    ("backend.routes.admin", None, ["admin"]),
    ("backend.routes.convert", None, ["convert"]),
    ("backend.routes.thumbnails", None, ["thumbnails"]),
    ("backend.routes.relax", "/api", ["relax"]),
    ("backend.routes.search", "/api/search", ["search"]),
    ("backend.routes.spectroscopy", "/api/spectroscopy", ["spectroscopy"]),
    ("backend.routes.energy", "/api/energy", ["energy"]),
    ("backend.routes.reaction", "/api/reaction", ["reaction"]),
    ("backend.routes.retrosynthesis", "/api/retrosynthesis", ["retrosynthesis"]),
    ("backend.routes.kab", "/api/kab", ["kab"]),
    ("backend.routes.quantum", "/api/quantum", ["quantum"]),
    ("backend.routes.collaboration", "/api/collaboration", ["collaboration"]),
    ("backend.routes.dashboard", "/api/dashboard", ["dashboard"]),
    ("backend.routes.screening", "/api/screening", ["screening"]),
    ("backend.routes.search_phase7", "/api/search", ["search-phase7"]),
    ("backend.routes.qm_md", "/api", ["qm-md"]),
    ("backend.api.search", "/api/search", ["search"]),
    ("backend.api.screening", "/api/screening", ["screening"]),
    ("backend.api.conformers", "/api/conformers", ["conformers"]),
    ("backend.api.molecule", "/api/molecule", ["molecule"]),
    ("backend.api.studio", "/api/studio", ["studio"]),
    ("backend.api.studio_v2", "/api/studio/v2", ["studio-v2"]),
]


def _should_skip_missing_dependency(exc: ModuleNotFoundError, module_path: str) -> bool:
    missing_name = exc.name or ""
    if not missing_name:
        return False
    if missing_name == module_path:
        return False
    return not missing_name.startswith("backend")


def _include_router(module_path: str, prefix: Optional[str], tags: list[str]) -> None:
    try:
        module = import_module(module_path)
    except ModuleNotFoundError as exc:
        if not _should_skip_missing_dependency(exc, module_path):
            raise
        logger.warning(
            "Skipping router %s because optional dependency '%s' is unavailable: %s",
            module_path,
            exc.name,
            exc,
        )
        return

    router = getattr(module, "router", None)
    if router is None:
        logger.warning("Skipping router %s because it does not expose 'router'", module_path)
        return

    include_kwargs = {"tags": tags}
    if prefix:
        include_kwargs["prefix"] = prefix
    app.include_router(router, **include_kwargs)


for module_path, prefix, tags in ROUTER_SPECS:
    _include_router(module_path, prefix, tags)



@app.get("/")
def root():
    return {
        "message": "MolForge Backend API",
        "version": settings.API_VERSION
    }


@app.get("/health")
def health():
    print("DEBUG: HEALTH CHECK HIT")
    return {"status": "healthy"}

