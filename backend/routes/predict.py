"""
Prediction route handlers
"""
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List
from backend.models.predictor import get_predictor, set_last_predictions, get_last_predictions
from backend.utils.featurizer import featurize_smiles, validate_smiles

router = APIRouter()


class PredictIn(BaseModel):
    smiles: str


class PredictOut(BaseModel):
    properties: dict


class BulkPredictIn(BaseModel):
    smiles_list: List[str]


class BulkPredictOut(BaseModel):
    results: List[dict]


@router.post("/", response_model=PredictOut)
def predict(payload: PredictIn):
    """
    Predict molecular properties from SMILES string
    """
    if not validate_smiles(payload.smiles):
        raise HTTPException(status_code=400, detail="Invalid SMILES string")
    
    features = featurize_smiles(payload.smiles)
    if features is None:
        raise HTTPException(status_code=400, detail="Failed to featurize SMILES")
    
    predictor = get_predictor()
    properties = predictor.predict(features)
    
    # Store last predictions for debugging
    set_last_predictions(properties)
    
    return PredictOut(properties=properties)


@router.get("/last")
def predict_last():
    """
    Get last predictions (for debugging)
    """
    last = get_last_predictions()
    if last is None:
        raise HTTPException(status_code=404, detail="No predictions available")
    return {"properties": last}


@router.post("/bulk", response_model=BulkPredictOut)
def predict_bulk(payload: BulkPredictIn):
    """
    Predict properties for multiple SMILES strings
    """
    predictor = get_predictor()
    results = []
    
    for smiles in payload.smiles_list:
        if not validate_smiles(smiles):
            results.append({
                "smiles": smiles,
                "error": "Invalid SMILES string",
                "properties": None
            })
            continue
        
        features = featurize_smiles(smiles)
        if features is None:
            results.append({
                "smiles": smiles,
                "error": "Failed to featurize SMILES",
                "properties": None
            })
            continue
        
        try:
            properties = predictor.predict(features)
            results.append({
                "smiles": smiles,
                "error": None,
                "properties": properties
            })
        except Exception as e:
            results.append({
                "smiles": smiles,
                "error": str(e),
                "properties": None
            })
    
    return BulkPredictOut(results=results)

