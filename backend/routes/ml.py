"""
ML Prediction API endpoints
"""
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Dict, Any, Optional
from chem.ml.predict import predict_properties
from chem.ml.features import extract_features
from chem.ml.explain import explain_prediction, explain_retrosynthesis_step, explain_reaction_mechanism

router = APIRouter()


class MoleculeRequest(BaseModel):
    atoms: List[Dict[str, Any]]
    bonds: List[Dict[str, Any]]


class ExplainRequest(BaseModel):
    molecule: MoleculeRequest
    property_name: Optional[str] = None


class PredictWithEnginesRequest(BaseModel):
    molecule: MoleculeRequest
    engine_outputs: Optional[Dict[str, Any]] = None
    use_engine_features: bool = False


class RetrosynthesisExplainRequest(BaseModel):
    pathway_step: Dict[str, Any]
    reaction: Optional[Dict[str, Any]] = None


class MechanismExplainRequest(BaseModel):
    mechanism_step: Dict[str, Any]
    step_index: int


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


@router.post("/predict/with-engines")
async def predict_ml_with_engines(req: PredictWithEnginesRequest):
    """Predict molecular properties using ML with engine-derived features"""
    try:
        molecule = {
            "atoms": req.molecule.atoms,
            "bonds": req.molecule.bonds
        }
        return predict_properties(
            molecule,
            engine_outputs=req.engine_outputs,
            use_engine_features=req.use_engine_features
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"ML prediction with engines failed: {str(e)}")


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


@router.post("/explain")
async def explain_ml(req: ExplainRequest):
    """Generate AI explanation for ML prediction"""
    try:
        molecule = {
            "atoms": req.molecule.atoms,
            "bonds": req.molecule.bonds
        }
        return explain_prediction(molecule, req.property_name)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Explanation generation failed: {str(e)}")


@router.post("/explain/retrosynthesis")
async def explain_retrosynthesis(req: RetrosynthesisExplainRequest):
    """Generate explanation for retrosynthesis step"""
    try:
        explanation = explain_retrosynthesis_step(req.pathway_step, req.reaction)
        return {
            "explanation": explanation,
            "step": req.pathway_step.get('step', 0)
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Retrosynthesis explanation failed: {str(e)}")


@router.post("/explain/mechanism")
async def explain_mechanism(req: MechanismExplainRequest):
    """Generate explanation for reaction mechanism step"""
    try:
        explanation = explain_reaction_mechanism(req.mechanism_step, req.step_index)
        return {
            "explanation": explanation,
            "step_index": req.step_index
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Mechanism explanation failed: {str(e)}")

