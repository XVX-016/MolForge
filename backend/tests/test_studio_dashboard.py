
import pytest
import uuid
from fastapi.testclient import TestClient
from backend.app import app

client = TestClient(app)

def test_studio_dashboard_placeholder():
    """Verify that the dashboard returns a default molecule for the zero-UUID placeholder"""
    zero_uuid = "00000000-0000-0000-0000-000000000000"
    response = client.post(
        "/api/molecule/dashboard",
        json={"baseline_version_id": zero_uuid}
    )
    
    assert response.status_code == 200
    data = response.json()
    assert data["baseline"]["version_id"] == zero_uuid
    assert data["baseline"]["smiles"] == "CC(=O)OC1=CC=CC=C1C(=O)O"
    assert "properties" in data["baseline"]
    assert "optimization_context" in data
    assert len(data["optimization_context"]["available_rules"]) > 0
