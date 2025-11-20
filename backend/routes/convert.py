"""
SMILES to Molfile conversion endpoint
Uses RDKit to generate 3D coordinates from SMILES
"""
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import Optional

router = APIRouter()


class SMILESIn(BaseModel):
    smiles: str


class MolfileOut(BaseModel):
    molfile: str


@router.post("/convert/smiles", response_model=MolfileOut)
def convert_smiles(payload: SMILESIn):
    """
    Convert SMILES string to 3D molfile (V2000 format)
    Uses RDKit to generate 3D coordinates and optimize geometry
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        raise HTTPException(
            status_code=503,
            detail="RDKit is not installed. Please install it to use this endpoint."
        )

    smiles = payload.smiles.strip()
    if not smiles:
        raise HTTPException(status_code=400, detail="SMILES string cannot be empty")

    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise HTTPException(status_code=400, detail="Invalid SMILES string")

        # Add hydrogens
        mol = Chem.AddHs(mol)

        # Generate 3D coordinates
        # EmbedMolecule returns 0 on success, -1 on failure
        embed_result = AllChem.EmbedMolecule(mol, randomSeed=42)
        if embed_result == -1:
            # Try with different method if first fails
            embed_result = AllChem.EmbedMolecule(mol, useRandomCoords=True)
            if embed_result == -1:
                raise HTTPException(
                    status_code=500,
                    detail="Failed to generate 3D coordinates for molecule"
                )

        # Optimize geometry using MMFF (Merck Molecular Force Field)
        try:
            AllChem.MMFFOptimizeMolecule(mol)
        except Exception:
            # If MMFF fails, try UFF (Universal Force Field) as fallback
            try:
                AllChem.UFFOptimizeMolecule(mol)
            except Exception:
                # If both fail, continue with unoptimized coordinates
                pass

        # Convert to molfile (V2000 format)
        molfile = Chem.MolToMolBlock(mol)

        return MolfileOut(molfile=molfile)

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Error converting SMILES to molfile: {str(e)}"
        )


