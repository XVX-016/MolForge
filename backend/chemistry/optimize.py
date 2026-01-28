from rdkit import Chem
from rdkit.Chem import FilterCatalog, QED, Descriptors
from typing import List, Dict, Any, Optional
from .rules import get_applicable_rules
import logging

logger = logging.getLogger(__name__)

# Initialize FilterCatalog with PAINS filters
params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
catalog = FilterCatalog.FilterCatalog(params)

class OptimizationIssue:
    def __init__(self, title: str, description: str, severity: str = "medium"):
        self.title = title
        self.description = description
        self.severity = severity # "high" (Red), "medium" (Amber), "info" (Blue/White)

    def to_dict(self):
        return {
            "title": self.title,
            "description": self.description,
            "severity": self.severity
        }

def analyze_structure(mol: Chem.Mol) -> List[Dict[str, Any]]:
    """
    Analyzes structure for medicinal chemistry issues (PAINS, Lipinski, etc.)
    """
    issues = []
    
    # 1. PAINS Detection
    matches = catalog.GetMatches(mol)
    for match in matches:
        issues.append(OptimizationIssue(
            title="PAINS Alert",
            description=f"Structural alert detected: {match.GetDescription()}",
            severity="high"
        ).to_dict())
        
    # 2. Molecular Weight Alert
    mw = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    if mw > 500:
        issues.append(OptimizationIssue(
            title="High Molecular Weight",
            description=f"MW is {mw:.1f}, exceeding Lipinski's limit (500). Consider truncating non-essential chains.",
            severity="medium"
        ).to_dict())
        
    # 3. LogP Alert
    logp = Chem.rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
    if logp > 5.0:
        issues.append(OptimizationIssue(
            title="High Lipophilicity",
            description=f"LogP is {logp:.1f}. High lipophilicity may reduce solubility and increase off-target binding.",
            severity="medium"
        ).to_dict())
    elif logp < 0:
         issues.append(OptimizationIssue(
            title="Low Lipophilicity",
            description=f"LogP is {logp:.1f}. Very polar molecules may have poor membrane permeability.",
            severity="info"
        ).to_dict())

    # 4. TPSA Alert
    tpsa = Chem.rdMolDescriptors.CalcTPSA(mol)
    if tpsa > 140:
        issues.append(OptimizationIssue(
            title="High TPSA",
            description=f"TPSA is {tpsa:.1f} Å². High polar surface area may limit blood-brain barrier penetration.",
            severity="medium"
        ).to_dict())

    # 5. Advanced Scoring (Phase 5B)
    qed_score = QED.qed(mol)
    if qed_score < 0.2:
        issues.append(OptimizationIssue(
            title="Low Drug-likeness (QED)",
            description=f"QED score is {qed_score:.2f}. Molecule lacks characteristics common to oral drugs.",
            severity="medium"
        ).to_dict())

    # 6. Radar Scores (Normalized 0.0 - 1.0)
    # We return these for the frontend multi-objective visualization
    radar = {
        "drug_likeness": round(qed_score, 2),
        "complexity": round(min(Descriptors.BertzCT(mol) / 1000.0, 1.0), 2),
        "lipophilicity": round(max(0, 1.0 - abs(Descriptors.MolLogP(mol) - 2.5) / 2.5), 2), # Peak at 2.5
        "solubility": round(max(0, 1.0 - Descriptors.TPSA(mol) / 140.0), 2),
        "size": round(max(0, 1.0 - Descriptors.MolWt(mol) / 500.0), 2)
    }

    return issues, radar

def suggest_optimizations(mol: Chem.Mol) -> List[Dict[str, Any]]:
    """
    Suggests bioisosteric modifications based on property profiles.
    Uses the formal ChemistryRule registry.
    """
    # For now, we return all applicable rules to give the AI context.
    # We could narrow this down based on the same heuristics as before.
    rules = get_applicable_rules(Chem.MolToSmiles(mol))
    return [rule.dict() for rule in rules]
