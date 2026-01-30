from typing import List, Optional, Literal, Union, Dict
from pydantic import BaseModel, Field

# --- Shared Components ---

class KeyObservation(BaseModel):
    metric: str = Field(..., description="The specific metric observed (e.g., 'Binding Energy', 'RMSD')")
    value: str = Field(..., description="The value of the metric")
    interpretation: str = Field(..., description="Scientific interpretation of this value")

class Risk(BaseModel):
    type: Literal["STABILITY", "BINDING", "MODEL_EXTRAPOLATION", "DATA_QUALITY"]
    severity: Literal["LOW", "MEDIUM", "HIGH"]
    description: str = Field(..., description="Details of the identified risk")

class SuggestedNextStep(BaseModel):
    action: Literal["RUN_MD", "REFINE_DOCKING", "BUILD_QSAR", "REJECT_COMPOUND", "FORCEFIELD_MINIMIZATION"]
    rationale: str = Field(..., description="Why this step is recommended based on current results")

class ComparativeContext(BaseModel):
    reference_range: str = Field(..., description="Standard range for this metric (e.g., '< -7.0 kcal/mol')")
    relative_position: str = Field(..., description="How this result compares (e.g., 'Superior to baseline')")

# --- Specific Extensions ---

class DockingAnalysis(BaseModel):
    binding_quality_score: float = Field(..., description="Binding affinity in kcal/mol")
    binding_quality_class: Literal["STRONG", "MODERATE", "WEAK"]
    pose_consistency: Literal["HIGH", "MEDIUM", "LOW"]
    interaction_summary: List[str] = Field(..., description="Key interactions (e.g., 'H-bond with ASP123')")

class MDAnalysis(BaseModel):
    stability_rmsd_trend: Literal["STABLE", "DRIFTING", "UNSTABLE"]
    equilibration_time_ns: float
    energy_convergence: Literal["CONVERGED", "NOT_CONVERGED"]
    warnings: List[str]

class QSARAnalysis(BaseModel):
    prediction_activity_type: str = Field(..., description="e.g. 'IC50'")
    prediction_value: float
    prediction_unit: str = Field("nM")
    applicability_domain: Literal["INSIDE", "BORDERLINE", "OUTSIDE"]
    model_confidence: Literal["HIGH", "MEDIUM", "LOW"]

# --- Main Insight Schema ---

class GeminiInsight(BaseModel):
    """
    Structured output schema for Gemini Scientific Assistant.
    This artifact is immutable once generated for a workflow node.
    """
    node_id: str
    analysis_type: Literal["DOCKING", "MD", "QSAR", "PHARMACOPHORE"]
    
    summary: str = Field(..., description="High-level executive summary of the findings")
    scientific_context: str = Field(..., description="Background context relevant to this specific analysis")
    
    key_findings: List[KeyObservation]
    risks: List[Risk]
    
    # Polymorphic analysis details
    specific_analysis: Union[DockingAnalysis, MDAnalysis, QSARAnalysis, None] = None
    
    comparative_context: Optional[ComparativeContext] = None
    suggested_next_steps: List[SuggestedNextStep]
    
    confidence_score: float = Field(..., ge=0.0, le=1.0, description="Model's self-reported confidence in this interpretation")
