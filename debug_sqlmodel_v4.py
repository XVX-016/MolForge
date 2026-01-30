from typing import Optional, Dict, Any
from sqlmodel import SQLModel, Field, create_engine
import uuid

class T1(SQLModel, table=True):
    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)

class T2(SQLModel, table=True):
    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    t1_id: uuid.UUID = Field(foreign_key="t1.id")

class T3(SQLModel, table=True):
    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    t1_id: uuid.UUID = Field(foreign_key="t1.id")
    t2_id: uuid.UUID = Field(foreign_key="t2.id")

print("Starting table creation...")
engine = create_engine("sqlite:///:memory:")
SQLModel.metadata.create_all(engine)
print("Success!")
