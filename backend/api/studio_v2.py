from fastapi import APIRouter, Depends, HTTPException, BackgroundTasks
from sqlmodel import Session, select
from typing import List, Optional
import uuid
from rdkit import Chem

from backend.core.dependencies import get_db
from backend.chemistry.models import (
    Experiment, MoleculeVersion, MoleculeIdentity, 
    WorkflowNode, WorkflowNodeResult, WorkflowNodeState
)
from backend.models.db.molecule import Molecule # Legacy model
from backend.services.studio_orchestrator import dispatch_node_execution
from backend.services.audit_service import log_user_action

router = APIRouter(tags=["Studio V2"])

@router.post("/experiment/from-library/{library_id}", response_model=Experiment)
def create_from_library(library_id: int, db: Session = Depends(get_db)):
    """Import a molecule from the legacy library and start an experiment"""
    # 1. Fetch legacy molecule
    legacy_mol = db.get(Molecule, library_id)
    if not legacy_mol:
        raise HTTPException(status_code=404, detail="Library molecule not found")
    
    if not legacy_mol.smiles:
         raise HTTPException(status_code=400, detail="Molecule has no SMILES")

    # 2. Compute Identity
    mol = Chem.MolFromSmiles(legacy_mol.smiles)
    if not mol:
        raise HTTPException(status_code=400, detail="Invalid SMILES")
    
    inchikey = Chem.MolToInchiKey(mol)
    canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    
    # 3. Find or Create Identity
    identity_query = select(MoleculeIdentity).where(MoleculeIdentity.inchikey == inchikey)
    identity = db.exec(identity_query).first()
    
    if not identity:
        identity = MoleculeIdentity(inchikey=inchikey)
        db.add(identity)
        db.commit()
        db.refresh(identity)
        
    # 4. Find or Create Version
    # Simplistic versioning: Check if ANY version exists, otherwise create v1
    version_query = select(MoleculeVersion).where(MoleculeVersion.molecule_id == identity.id).order_by(MoleculeVersion.version_index.desc())
    latest_version = db.exec(version_query).first()
    
    if not latest_version:
        latest_version = MoleculeVersion(
            molecule_id=identity.id,
            version_index=1,
            canonical_smiles=canonical_smiles,
            state="VERSIONED",
            properties={"imported_from_library": library_id}
        )
        db.add(latest_version)
        db.commit()
        db.refresh(latest_version)
        
    # 5. Delegate to create_experiment
    return create_experiment(molecule_version_id=latest_version.id, forcefield="GAFF2", db=db)

@router.post("/experiment", response_model=Experiment)
def create_experiment(
    molecule_version_id: uuid.UUID,
    forcefield: str = "GAFF2",
    db: Session = Depends(get_db)
):
    """Create a new drug discovery experiment snapshot"""
    version = db.get(MoleculeVersion, molecule_version_id)
    if not version:
        raise HTTPException(status_code=404, detail="Molecule version not found")
    
    experiment = Experiment(
        molecule_id=version.molecule_id,
        molecule_version_id=molecule_version_id,
        forcefield=forcefield
    )
    db.add(experiment)
    db.commit()
    db.refresh(experiment)
    
    # Automatically create the Root/Baseline Node
    # This ensures the timeline always has a starting point (The "Ghost Node" becomes real)
    root_node = WorkflowNode(
        experiment_id=experiment.id,
        node_type="BASELINE",
        parameters_hash="ROOT",  # Root is unique by definition
        input_params={"description": "Experiment Root", "version": version.version_index}
    )
    db.add(root_node)
    
    # Also create the initial Result for the root node so it shows as COMPLETED
    root_result = WorkflowNodeResult(
        node_id=root_node.id,
        status=WorkflowNodeState.COMPLETED,
        output_data={"message": "Experiment initialized"}
    )
    db.add(root_result)
    
    # Audit Log
    log_user_action(
        db, 
        action_type="CREATE_EXPERIMENT", 
        entity_type="Experiment", 
        entity_id=str(experiment.id),
        details={
            "molecule_id": str(version.molecule_id),
            "version_index": version.version_index,
            "forcefield": forcefield
        }
    )
    
    db.commit()
    
    return experiment

@router.get("/experiment/{experiment_id}", response_model=Experiment)
def get_experiment(experiment_id: uuid.UUID, db: Session = Depends(get_db)):
    experiment = db.get(Experiment, experiment_id)
    if not experiment:
        raise HTTPException(status_code=404, detail="Experiment not found")
    return experiment

@router.post("/experiment/{experiment_id}/node", response_model=WorkflowNode)
def add_workflow_node(
    experiment_id: uuid.UUID,
    node_type: str,
    input_params: dict,
    parent_node_id: Optional[uuid.UUID] = None,
    db: Session = Depends(get_db)
):
    """Add an intent node to the workflow timeline"""
    experiment = db.get(Experiment, experiment_id)
    if not experiment:
        raise HTTPException(status_code=404, detail="Experiment not found")
    
    # Placeholder for parameter hashing
    param_hash = str(hash(frozenset(input_params.items())))

    # Strict Governance: Ensure parent node belongs to the same experiment (DAG integrity)
    if parent_node_id:
        parent_node = db.get(WorkflowNode, parent_node_id)
        if not parent_node:
            raise HTTPException(status_code=404, detail="Parent node not found")
        if parent_node.experiment_id != experiment_id:
            raise HTTPException(status_code=400, detail="Parent node belongs to a different experiment")

    node = WorkflowNode(
        experiment_id=experiment_id,
        node_type=node_type,
        parameters_hash=param_hash,
        input_params=input_params,
        parent_node_id=parent_node_id
    )
    db.add(node)

    log_user_action(
        db,
        action_type="ADD_NODE",
        entity_type="WorkflowNode",
        entity_id=str(node.id),
        details={
            "node_type": node_type,
            "parent_node_id": str(parent_node_id) if parent_node_id else None,
            "experiment_id": str(experiment_id)
        }
    )

    db.commit()
    db.refresh(node)
    return node

@router.post("/node/{node_id}/run", response_model=WorkflowNodeResult)
def run_node(
    node_id: uuid.UUID,
    background_tasks: BackgroundTasks,
    db: Session = Depends(get_db)
):
    """Trigger the execution of a scientific node"""
    node = db.get(WorkflowNode, node_id)
    if not node:
        raise HTTPException(status_code=404, detail="Workflow node not found")
    
    # Create the result entry (The execution record)
    result = WorkflowNodeResult(
        node_id=node_id,
        status=WorkflowNodeState.QUEUED
    )
    db.add(result)

    log_user_action(
        db,
        action_type="RUN_NODE",
        entity_type="WorkflowNode",
        entity_id=str(node_id),
        details={"result_id": str(result.id)}
    )

    db.commit()
    db.refresh(result)
    
    # Dispatch background task
    background_tasks.add_task(dispatch_node_execution, node_id, result.id)
    
    return result

@router.get("/experiment/{experiment_id}/workflow", response_model=List[WorkflowNode])
def get_workflow_timeline(experiment_id: uuid.UUID, db: Session = Depends(get_db)):
    """Get the full sequence of analysis nodes for an experiment"""
    statement = select(WorkflowNode).where(WorkflowNode.experiment_id == experiment_id)
    results = db.exec(statement).all()
    return results
