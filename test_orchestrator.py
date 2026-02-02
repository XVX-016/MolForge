import asyncio
import uuid
import sys
import os
from sqlmodel import Session, select
from backend.db import engine, init_db
from backend.chemistry.models import (
    MoleculeIdentity, MoleculeVersion, Experiment, 
    WorkflowNode, WorkflowNodeResult, WorkflowNodeState
)
from backend.services.studio_orchestrator import dispatch_node_execution

async def test_workflow():
    print("Starting integration test...")
    init_db()
    
    with Session(engine) as db:
        # 1. Setup Sample Data
        molecule = MoleculeIdentity(inchikey="TESTKEY" + str(uuid.uuid4())[:8])
        db.add(molecule)
        db.commit()
        db.refresh(molecule)
        
        version = MoleculeVersion(
            molecule_id=molecule.id,
            canonical_smiles="CC",
            version_index=1
        )
        db.add(version)
        db.commit()
        db.refresh(version)
        
        # 2. Create Experiment
        experiment = Experiment(
            molecule_id=molecule.id,
            molecule_version_id=version.id
        )
        db.add(experiment)
        db.commit()
        db.refresh(experiment)
        print(f"Created Experiment: {experiment.id}")
        
        # 3. Add Node
        node = WorkflowNode(
            experiment_id=experiment.id,
            node_type="docking",
            parameters_hash="test_hash",
            input_params={"target": "protein.pdb"}
        )
        db.add(node)
        db.commit()
        db.refresh(node)
        print(f"Added Node: {node.id}")
        
        # 4. Initialize Result
        result = WorkflowNodeResult(
            node_id=node.id,
            status=WorkflowNodeState.QUEUED
        )
        db.add(result)
        db.commit()
        db.refresh(result)
        
        # 5. Run Orchestrator
        print("Dispatching orchestrator...")
        await dispatch_node_execution(node.id, result.id)
        
        # 6. Verify Completion
        db.refresh(result)
        print(f"Result Status: {result.status}")
        print(f"Output Data: {result.output_data}")
        
        assert result.status == WorkflowNodeState.COMPLETED
        assert "binding_affinity" in result.output_data
        print("Integration Test SUCCESS!")

if __name__ == "__main__":
    asyncio.run(test_workflow())
