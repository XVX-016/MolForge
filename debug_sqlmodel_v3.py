from typing import Optional, Dict, Any, List
from datetime import datetime
from sqlmodel import SQLModel, Field, JSON, Column, create_engine
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
    __tablename__ = "molecules"
    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    inchikey: str = Field(unique=True, index=True)
    created_at: datetime = Field(default_factory=datetime.utcnow)
    user_id: Optional[uuid.UUID] = Field(default=None, index=True)

class MoleculeVersion(SQLModel, table=True):
    __tablename__ = "molecule_versions"
    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    molecule_id: uuid.UUID = Field(foreign_key="molecules.id", index=True)
    version_index: int = Field(default=1)
    parent_version_id: Optional[uuid.UUID] = Field(default=None)
    canonical_smiles: str = Field(index=True)
    json_graph: Dict[str, Any] = Field(default_factory=dict, sa_column=Column(JSON))
    properties: Dict[str, Any] = Field(default_factory=dict, sa_column=Column(JSON))
    state: str = Field(default="VERSIONED")
    additional_metadata: Dict[str, Any] = Field(default_factory=dict, sa_column=Column(JSON))
    created_at: datetime = Field(default_factory=datetime.utcnow)
    user_id: Optional[uuid.UUID] = Field(default=None, index=True)

class Experiment(SQLModel, table=True):
    __tablename__ = "experiments"
    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    molecule_id: uuid.UUID = Field(index=True) # Redundant foreign key removed to avoid mapper ambiguity
    molecule_version_id: uuid.UUID = Field(foreign_key="molecule_versions.id")
    forcefield: str = Field(default="GAFF2")
    created_at: datetime = Field(default_factory=datetime.utcnow)
    updated_at: datetime = Field(default_factory=datetime.utcnow)
    metadata: Dict[str, Any] = Field(default_factory=dict, sa_column=Column(JSON))

class WorkflowNode(SQLModel, table=True):
    __tablename__ = "workflow_nodes"
    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    experiment_id: uuid.UUID = Field(foreign_key="experiments.id", index=True)
    node_type: str = Field(index=True)
    parameters_hash: str = Field(index=True)
    input_params: Dict[str, Any] = Field(default_factory=dict, sa_column=Column(JSON))
    created_at: datetime = Field(default_factory=datetime.utcnow)

class WorkflowNodeResult(SQLModel, table=True):
    __tablename__ = "workflow_node_results"
    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    node_id: uuid.UUID = Field(foreign_key="workflow_nodes.id", index=True)
    status: WorkflowNodeState = Field(default=WorkflowNodeState.CREATED)
    stdout_logs: Optional[str] = None
    output_data: Dict[str, Any] = Field(default_factory=dict, sa_column=Column(JSON))
    artifacts_checksum: Optional[str] = None
    created_at: datetime = Field(default_factory=datetime.utcnow)

print("Starting table creation...")
engine = create_engine("sqlite:///:memory:")
SQLModel.metadata.create_all(engine)
print("Success!")
