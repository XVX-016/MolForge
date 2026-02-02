from typing import Optional, Dict, Any, List
from datetime import datetime
from sqlmodel import SQLModel, Field, JSON, Column
import uuid
from enum import Enum

class WorkflowNodeState(str, Enum):
    CREATED = "created"
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    INVALIDATED = "invalidated"

class MoleculeIdentity(SQLModel, table=True):
    """
    The conceptual identity of a molecule. 
    Unique per InChIKey.
    """
    __tablename__ = "molecules"
    
    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    inchikey: str = Field(unique=True, index=True)
    created_at: datetime = Field(default_factory=datetime.utcnow)
    user_id: Optional[uuid.UUID] = Field(default=None, index=True)

class MoleculeVersion(SQLModel, table=True):
    """
    An immutable snapshot of a molecule at a specific point in time.
    """
    __tablename__ = "molecule_versions"
    
    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    molecule_id: uuid.UUID = Field(foreign_key="molecules.id", index=True)
    version_index: int = Field(default=1)
    parent_version_id: Optional[uuid.UUID] = Field(default=None, foreign_key="molecule_versions.id")
    
    # Snapshotted Data
    canonical_smiles: str = Field(index=True)
    json_graph: Dict[str, Any] = Field(default_factory=dict, sa_column=Column(JSON))
    properties: Dict[str, Any] = Field(default_factory=dict, sa_column=Column(JSON))
    
    # Lifecycle State
    state: str = Field(default="VERSIONED") # VERSIONED | OPTIMIZED
    
    # Additional metadata
    additional_metadata: Dict[str, Any] = Field(default_factory=dict, sa_column=Column(JSON))
    created_at: datetime = Field(default_factory=datetime.utcnow)
    user_id: Optional[uuid.UUID] = Field(default=None, index=True)

class Experiment(SQLModel, table=True):
    """
    A root container for a drug discovery workflow session.
    """
    __tablename__ = "experiments"
    
    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    molecule_id: uuid.UUID = Field(foreign_key="molecules.id", index=True)
    molecule_version_id: uuid.UUID = Field(foreign_key="molecule_versions.id")
    
    # Snapshot parameters
    forcefield: str = Field(default="GAFF2")
    
    created_at: datetime = Field(default_factory=datetime.utcnow)
    updated_at: datetime = Field(default_factory=datetime.utcnow)
    additional_metadata: Dict[str, Any] = Field(default_factory=dict, sa_column=Column(JSON))

class WorkflowNode(SQLModel, table=True):
    """
    The intent to perform a scientific calculation.
    """
    __tablename__ = "workflow_nodes"
    
    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    experiment_id: uuid.UUID = Field(foreign_key="experiments.id", index=True)
    node_type: str = Field(index=True) # docking, md, etc.
    parent_node_id: Optional[uuid.UUID] = Field(default=None, foreign_key="workflow_nodes.id", index=True)
    
    # Parameters hash for deduplication/rerun logic
    parameters_hash: str = Field(index=True)
    input_params: Dict[str, Any] = Field(default_factory=dict, sa_column=Column(JSON))
    
    created_at: datetime = Field(default_factory=datetime.utcnow)

class WorkflowNodeResult(SQLModel, table=True):
    """
    The execution immutable artifact of a WorkflowNode.
    """
    __tablename__ = "workflow_node_results"
    
    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    node_id: uuid.UUID = Field(foreign_key="workflow_nodes.id", index=True)
    
    status: WorkflowNodeState = Field(default=WorkflowNodeState.CREATED)
    stdout_logs: Optional[str] = None
    
    # Artifacts and scientific data
    output_data: Dict[str, Any] = Field(default_factory=dict, sa_column=Column(JSON))
    artifacts_checksum: Optional[str] = None # SHA-256 of generated files
    
    created_at: datetime = Field(default_factory=datetime.utcnow)
