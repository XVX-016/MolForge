"""
Quantum Properties API endpoints
"""
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Dict, Any
from chem.quantum.homo_lumo import calculate_homo_lumo
from chem.quantum.esp import calculate_esp

router = APIRouter()


class MoleculeRequest(BaseModel):
    atoms: List[Dict[str, Any]]
    bonds: List[Dict[str, Any]]


@router.post("/calculate")
async def calculate_quantum(req: MoleculeRequest):
    """Calculate HOMO/LUMO and ESP"""
    try:
        molecule = {
            "atoms": req.atoms,
            "bonds": req.bonds
        }
        homo_lumo = calculate_homo_lumo(molecule)
        esp = calculate_esp(molecule)
        return {
            "homo_lumo": homo_lumo,
            "esp": esp
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Quantum calculation failed: {str(e)}")


@router.post("/homo-lumo")
async def calculate_homo_lumo_route(req: MoleculeRequest):
    """Calculate HOMO/LUMO orbitals"""
    try:
        molecule = {
            "atoms": req.atoms,
            "bonds": req.bonds
        }
        return calculate_homo_lumo(molecule)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"HOMO/LUMO calculation failed: {str(e)}")


@router.post("/esp")
async def calculate_esp_route(req: MoleculeRequest):
    """Calculate electrostatic potential"""
    try:
        molecule = {
            "atoms": req.atoms,
            "bonds": req.bonds
        }
        return calculate_esp(molecule)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"ESP calculation failed: {str(e)}")

