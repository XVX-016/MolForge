"""
Test prediction endpoints
"""
import pytest
from fastapi.testclient import TestClient
from backend.app import app

client = TestClient(app)


def test_predict_valid_smiles():
    """Test prediction with valid SMILES"""
    resp = client.post("/predict/", json={"smiles": "CCO"})
    assert resp.status_code == 200
    assert "properties" in resp.json()
    
    props = resp.json()["properties"]
    assert "stability" in props
    assert "toxicity" in props
    assert "solubility" in props
    assert "bioavailability" in props
    assert "novelty" in props
    
    # Check that values are floats
    assert isinstance(props["stability"], float)
    assert isinstance(props["toxicity"], float)


def test_predict_invalid_smiles():
    """Test prediction with invalid SMILES"""
    resp = client.post("/predict/", json={"smiles": "INVALID"})
    assert resp.status_code == 400
    assert "detail" in resp.json()


def test_predict_last():
    """Test getting last predictions"""
    # First make a prediction
    client.post("/predict/", json={"smiles": "CCO"})
    
    # Then get last predictions
    resp = client.get("/predict/last")
    assert resp.status_code == 200
    assert "properties" in resp.json()


def test_predict_bulk():
    """Test bulk prediction"""
    resp = client.post("/predict/bulk", json={
        "smiles_list": ["CCO", "CC", "C"]
    })
    assert resp.status_code == 200
    assert "results" in resp.json()
    assert len(resp.json()["results"]) == 3


def test_predict_fast():
    """Test fast prediction endpoint"""
    resp = client.post("/predict-fast", json={"smiles": "CCO"})
    # Should work or fallback to regular predict
    assert resp.status_code in [200, 400]  # 400 if ONNX not available
    if resp.status_code == 200:
        assert "properties" in resp.json()

