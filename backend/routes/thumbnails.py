"""
Thumbnail generation endpoints
"""
from fastapi import APIRouter, HTTPException
from fastapi.responses import Response
from pydantic import BaseModel
from typing import Optional, Tuple
from backend.services.thumbnail_service import ThumbnailService

router = APIRouter()


class ThumbnailRequest(BaseModel):
    """Request model for thumbnail generation"""
    molfile: Optional[str] = None
    smiles: Optional[str] = None
    size: Optional[Tuple[int, int]] = (600, 600)
    format: Optional[str] = "PNG"


class ThumbnailBase64Request(BaseModel):
    """Request model for base64 thumbnail generation"""
    molfile: Optional[str] = None
    smiles: Optional[str] = None
    size: Optional[Tuple[int, int]] = (600, 600)


@router.post("/thumbnails/generate")
async def generate_thumbnail(request: ThumbnailRequest):
    """
    Generate a 2D molecule thumbnail image
    
    Accepts either molfile or SMILES string.
    Returns PNG image bytes.
    """
    try:
        img_bytes = ThumbnailService.generate_2d_thumbnail(
            molfile=request.molfile,
            smiles=request.smiles,
            size=request.size,
            format=request.format
        )
        
        return Response(
            content=img_bytes,
            media_type="image/png",
            headers={"Content-Disposition": "inline; filename=thumbnail.png"}
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error generating thumbnail: {str(e)}")


@router.post("/thumbnails/generate-base64")
async def generate_thumbnail_base64(request: ThumbnailBase64Request):
    """
    Generate a 2D molecule thumbnail and return as base64 string
    
    Accepts either molfile or SMILES string.
    Returns base64-encoded image (data URI format).
    """
    try:
        thumbnail_b64 = ThumbnailService.generate_2d_thumbnail_base64(
            molfile=request.molfile,
            smiles=request.smiles,
            size=request.size
        )
        
        return {"thumbnail_b64": thumbnail_b64}
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error generating thumbnail: {str(e)}")


@router.post("/thumbnails/generate-3d")
async def generate_3d_molfile(request: ThumbnailBase64Request):
    """
    Generate optimized 3D molfile from SMILES or existing molfile
    
    Accepts either molfile or SMILES string.
    Returns optimized 3D molfile.
    """
    try:
        molfile_3d = ThumbnailService.generate_3d_molfile(
            molfile=request.molfile,
            smiles=request.smiles
        )
        
        return {"molfile": molfile_3d}
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error generating 3D structure: {str(e)}")

