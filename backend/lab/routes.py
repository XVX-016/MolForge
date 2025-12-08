from fastapi import APIRouter
from pydantic import BaseModel
from typing import List, Tuple

router = APIRouter()

class AtomIn(BaseModel):
    id: str
    element: str
    position: Tuple[float,float,float]

class PredictRequest(BaseModel):
    atoms: List[AtomIn]

class BondOut(BaseModel):
    id: str
    a: str
    b: str
    order: int

class PredictResponse(BaseModel):
    bonds: List[BondOut]

# import predictor
from backend.chem.ml.bond_predictor import BondPredictor
predictor = BondPredictor()

@router.post('/predict-bonds', response_model=PredictResponse)
async def predict_bonds(req: PredictRequest):
    atoms = [ (a.id, a.element, a.position) for a in req.atoms ]
    bonds = predictor.predict(atoms)
    return PredictResponse(bonds=bonds)

@router.post('/parse-molecule')
async def parse_molecule():
    # TODO: Implement parse-molecule
    return {"status": "not implemented"}

@router.post('/save-session')
async def save_session():
    # TODO: Implement save-session
    return {"status": "not implemented"}

@router.post('/run-ml-prediction')
async def run_ml_prediction():
    # TODO: Implement run-ml-prediction
    return {"status": "not implemented"}
