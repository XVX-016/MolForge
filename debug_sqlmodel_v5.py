from typing import Optional, Dict, Any, List
from datetime import datetime
from sqlmodel import SQLModel, Field, JSON, Column, create_engine
import uuid
from enum import Enum

class WorkflowNodeState(str, Enum):
    CREATED = "created"
    COMPLETED = "completed"

class MoleculeIdentity(SQLModel, table=True):
    __tablename__ = "molecules"
    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)

class MoleculeVersion(SQLModel, table=True):
    __tablename__ = "molecule_versions"
    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    molecule_id: uuid.UUID = Field(foreign_key="molecules.id", index=True)
    json_graph: Dict[str, Any] = Field(default_factory=dict, sa_column=Column(JSON))

class Experiment(SQLModel, table=True):
    __tablename__ = "experiments"
    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    molecule_id: uuid.UUID = Field(foreign_key="molecules.id", index=True)
    molecule_version_id: uuid.UUID = Field(foreign_key="molecule_versions.id")

print("Starting table creation...")
engine = create_engine("sqlite:///:memory:")
SQLModel.metadata.create_all(engine)
print("Success!")
