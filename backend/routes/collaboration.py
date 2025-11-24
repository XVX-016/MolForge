"""
Collaboration & Cloud Storage API endpoints
"""
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Dict, Any, Optional
from chem.collaboration.cloud import save_molecule_cloud, load_molecule_cloud, list_user_molecules
from chem.collaboration.versioning import create_version, get_version_history, fork_molecule

router = APIRouter()


class MoleculeRequest(BaseModel):
    atoms: List[Dict[str, Any]]
    bonds: List[Dict[str, Any]]


class SaveMoleculeRequest(BaseModel):
    molecule: MoleculeRequest
    user_id: str
    molecule_name: str
    metadata: Optional[Dict[str, Any]] = None


class ForkRequest(BaseModel):
    source_molecule_id: str
    new_owner_id: str
    fork_name: Optional[str] = None


@router.post("/save")
async def save_molecule(req: SaveMoleculeRequest):
    """Save molecule to cloud storage"""
    try:
        molecule = {
            "atoms": req.molecule.atoms,
            "bonds": req.molecule.bonds
        }
        return save_molecule_cloud(molecule, req.user_id, req.molecule_name, req.metadata)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Save failed: {str(e)}")


@router.get("/load/{molecule_id}")
async def load_molecule(molecule_id: str, user_id: Optional[str] = None):
    """Load molecule from cloud storage"""
    try:
        result = load_molecule_cloud(molecule_id, user_id)
        if result is None:
            raise HTTPException(status_code=404, detail="Molecule not found")
        return result
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Load failed: {str(e)}")


@router.get("/list/{user_id}")
async def list_molecules(user_id: str, limit: int = 50):
    """List user's molecules"""
    try:
        return list_user_molecules(user_id, limit)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"List failed: {str(e)}")


@router.post("/version")
async def create_version_route(molecule_id: str, molecule: MoleculeRequest, editor_id: str, description: Optional[str] = None):
    """Create new molecule version"""
    try:
        mol_dict = {
            "atoms": molecule.atoms,
            "bonds": molecule.bonds
        }
        return create_version(molecule_id, mol_dict, editor_id, description)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Version creation failed: {str(e)}")


@router.get("/history/{molecule_id}")
async def get_history(molecule_id: str):
    """Get molecule version history"""
    try:
        return get_version_history(molecule_id)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"History retrieval failed: {str(e)}")


@router.post("/fork")
async def fork_molecule_route(req: ForkRequest):
    """Fork a molecule"""
    try:
        return fork_molecule(req.source_molecule_id, req.new_owner_id, req.fork_name)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Fork failed: {str(e)}")

