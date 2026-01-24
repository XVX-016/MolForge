"""
Studio API Endpoints
"""
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import Dict, Any, Optional
import logging
from backend.services.studio_service import get_studio_service

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/studio", tags=["studio"])

class StudioCommandRequest(BaseModel):
    prompt: str
    molecule_context: Dict[str, Any]
    mode: str
    analysis_context: Optional[Dict[str, Any]] = None

@router.post("/command")
async def process_studio_command(request: StudioCommandRequest):
    """
    Process a Studio AI command.
    """
    # DEBUG: PROVE THIS CODE IS RUNNING
    print(f"DEBUG: HIT STUDIO API with prompt: {request.prompt}")
    
    try:
        service = get_studio_service()
        action = await service.process_command(
            request.prompt, 
            request.molecule_context, 
            request.mode,
            request.analysis_context
        )
        return action
    except Exception as e:
        logger.error(f"Studio API Error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))
