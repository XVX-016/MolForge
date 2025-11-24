"""
ML Prediction API endpoints
"""
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Dict, Any
from chem.ml.predict import predict_properties
from chem.ml.features import extract_features

router = APIRouter()


class MoleculeRequest(BaseModel):
    atoms: List[Dict[str, Any]]
    bonds: List[Dict[str, Any]]


@router.post("/predict")
async def predict_ml(req: MoleculeRequest):
    """Predict molecular properties using ML"""
    try:
        molecule = {
            "atoms": req.atoms,
            "bonds": req.bonds
        }
        return predict_properties(molecule)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"ML prediction failed: {str(e)}")


@router.post("/features")
async def extract_features_route(req: MoleculeRequest):
    """Extract molecular features"""
    try:
        molecule = {
            "atoms": req.atoms,
            "bonds": req.bonds
        }
        return extract_features(molecule)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Feature extraction failed: {str(e)}")

