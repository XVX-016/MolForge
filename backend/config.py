
import os
from dotenv import load_dotenv
from typing import Optional, List
from pathlib import Path
env_path = Path(__file__).resolve().parent / ".env"
load_dotenv(dotenv_path=env_path, override=True)


class Settings:
    DATABASE_URL: str = os.getenv(
        "DATABASE_URL",
        "sqlite:///./biosynth.db"
    )
    API_TITLE: str = "MolForge Backend API"
    API_VERSION: str = "0.1.0"
    API_DESCRIPTION: str = "Backend API for MolForge molecular design platform"
    _cors_origins_env = os.getenv("CORS_ORIGINS", "")
    _is_dev = os.getenv("ENVIRONMENT", "development").lower() in ["development", "dev"]
    
    if _cors_origins_env:
        CORS_ORIGINS: List[str] = [origin.strip() for origin in _cors_origins_env.split(",") if origin.strip()]
        CORS_ALLOW_ORIGIN_REGEX: Optional[str] = None
    elif _is_dev:
        CORS_ORIGINS: List[str] = [
            "http://localhost:5173",
            "http://localhost:5174",
            "http://localhost:3000",
            "http://localhost:5175",
            "http://localhost:8080",
            "http://127.0.0.1:5173",
            "http://127.0.0.1:5174",
            "http://127.0.0.1:3000",
            "http://127.0.0.1:5175",
            "http://127.0.0.1:8080",
        ]
        CORS_ALLOW_ORIGIN_REGEX: Optional[str] = r"^https?://(localhost|127\.0\.0\.1)(:\d+)?$"
    else:
        CORS_ORIGINS: List[str] = [
            "http://localhost:5173",
            "http://localhost:5174",
            "http://127.0.0.1:5173",
            "http://127.0.0.1:5174",
        ]
        CORS_ALLOW_ORIGIN_REGEX: Optional[str] = None
    MODEL_WEIGHTS_PATH: str = os.getenv(
        "MODEL_WEIGHTS_PATH",
        "backend/weights/property_predictor.pt"
    )
    
    ONNX_MODEL_PATH: str = os.getenv(
        "ONNX_MODEL_PATH",
        "backend/weights/property_predictor.onnx"
    )
    SECRET_KEY: Optional[str] = os.getenv("SECRET_KEY")
    ALGORITHM: str = "HS256"
    ACCESS_TOKEN_EXPIRE_MINUTES: int = 30
    LOG_LEVEL: str = os.getenv("LOG_LEVEL", "INFO")
    
    # AI Keys
    GEMINI_KEY: Optional[str] = os.getenv("GEMINI_KEY")
    GOOGLE_API_KEY: Optional[str] = os.getenv("GOOGLE_API_KEY")

settings = Settings()

