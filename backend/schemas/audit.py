
from pydantic import BaseModel, Field
from typing import List, Dict, Any, Optional
from datetime import datetime
import uuid

class AIAuditLog(BaseModel):
    """
    Regulatory-grade audit log for AI-assisted molecular design.
    Ensures 21 CFR Part 11 compliance by tracking all non-authoritative AI inputs.
    """
    id: uuid.UUID = Field(default_factory=uuid.uuid4)
    timestamp: datetime = Field(default_factory=datetime.utcnow)
    
    # AI Environment
    model_provider: str = Field(..., description="e.g., 'google'")
    model_name: str = Field(..., description="e.g., 'gemini-2.0-flash-exp'")
    system_prompt_hash: str = Field(..., description="SHA-256 hash of the system instructions at time of execution")
    
    # Input/Output
    input_messages_hash: str = Field(..., description="SHA-256 hash of the full prompt sent to AI")
    output_json: Dict[str, Any] = Field(..., description="Raw JSON command received from AI")
    
    # Chemistry Context
    baseline_version_id: uuid.UUID = Field(..., description="The molecule state Gemini was planning against")
    rule_ids_used: List[str] = Field(..., description="Which deterministic rules were extracted from the AI response")

    # Accountability
    human_signoff_id: Optional[str] = Field(None, description="ID of the chemist who approved this proposal")
    signoff_timestamp: Optional[datetime] = Field(None)
    signoff_reason: Optional[str] = Field(None)

    class Config:
        json_encoders = {
            datetime: lambda v: v.isoformat(),
            uuid.UUID: lambda v: str(v)
        }
