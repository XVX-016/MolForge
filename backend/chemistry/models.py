from typing import Optional, Dict, Any, List
from datetime import datetime
from sqlmodel import SQLModel, Field, JSON, Column
import uuid

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
