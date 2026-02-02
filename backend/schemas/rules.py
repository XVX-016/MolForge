
from pydantic import BaseModel, Field
from typing import List, Dict, Any, Optional
from enum import Enum

class RuleOperation(str, Enum):
    REPLACE = "REPLACE"
    ADD = "ADD"
    REMOVE = "REMOVE"
    PROTECT = "PROTECT"

class PropertyImpact(BaseModel):
    logp_delta: float = Field(default=0.0, description="Estimated change in LogP")
    tpsa_delta: float = Field(default=0.0, description="Estimated change in TPSA (Å²)")
    mw_delta: float = Field(default=0.0, description="Estimated change in Molecular Weight")
    qed_delta: float = Field(default=0.0, description="Estimated change in QED score")

class RuleMetadata(BaseModel):
    rationale: str = Field(..., description="Scientific motivation for this rule")
    typical_use: str = Field(..., description="When a chemist would typically apply this")
    tags: List[str] = Field(default_factory=list, description="Categorization tags (e.g., 'solubility', 'metabolic-stability')")

class ChemistryRule(BaseModel):
    id: str = Field(..., description="Unique identifier for the rule (e.g., 'RING_CLOSURE')")
    title: str = Field(..., description="Human-readable title")
    description: str = Field(..., description="Detailed description of what the rule does")
    operation: RuleOperation = Field(..., description="Type of chemical operation")
    
    # RDKit execution parameters
    smarts_pattern: str = Field(..., description="SMARTS pattern to match the target site")
    replacement_smiles: Optional[str] = Field(None, description="SMILES to use for REPLACE or ADD operations")
    
    # AI Reasoning context
    impact: PropertyImpact = Field(default_factory=PropertyImpact)
    metadata: RuleMetadata = Field(...)

    class Config:
        schema_extra = {
            "example": {
                "id": "opt-pyridine",
                "title": "Phenyl to Pyridine",
                "description": "Replace a phenyl ring with a pyridine ring to reduce LogP and increase solubility.",
                "operation": "REPLACE",
                "smarts_pattern": "c1ccccc1",
                "replacement_smiles": "c1ccncc1",
                "impact": {
                    "logp_delta": -0.8,
                    "tpsa_delta": 12.9,
                    "mw_delta": 1.0
                },
                "metadata": {
                    "rationale": "Inserting a nitrogen atom reduces lipophilicity and increases polar character.",
                    "typical_use": "When LogP is too high and compound has poor solubility.",
                    "tags": ["solubility", "logp", "bioisostere"]
                }
            }
        }

class RuleSet(BaseModel):
    version: str = "1.0.0"
    rules: List[ChemistryRule]

class ProposalStep(BaseModel):
    rule_id: str = Field(..., description="ID of the rule to apply")
    reason: str = Field(..., description="Justification for this specific step")
    priority: int = Field(default=1, description="Execution order priority")

class ProposalGraph(BaseModel):
    """
    Formalized intent graph proposed by Gemini.
    """
    baseline_smiles: str = Field(..., description="The starting material for this plan")
    steps: List[ProposalStep] = Field(..., max_items=3, description="Sequential optimization steps")
    rationale: str = Field(..., description="Overall strategy rationale")
