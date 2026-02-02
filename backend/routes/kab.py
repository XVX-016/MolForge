"""
KAB (Knowledge, Analysis & Binding) API endpoints
"""
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Dict, Any, Optional
from backend.chem.kab.binding import predict_binding_sites, score_sites
from backend.chem.kab.analysis import compute_similarity, substructure_alerts, analyze_molecule
from backend.chem.kab.knowledge import query_knowledge_rules, apply_ml_rules

router = APIRouter()


class MoleculeRequest(BaseModel):
    atoms: List[Dict[str, Any]]
    bonds: List[Dict[str, Any]]


class SimilarityRequest(BaseModel):
    molecule1: MoleculeRequest
    molecule2: MoleculeRequest


class ExplainRequest(BaseModel):
    molecule: MoleculeRequest
    site_id: Optional[int] = None
    alert_id: Optional[int] = None
    rule_id: Optional[str] = None


@router.post("/binding")
async def predict_binding_route(req: MoleculeRequest):
    """Predict binding sites for molecule"""
    try:
        molecule = {
            "atoms": req.atoms,
            "bonds": req.bonds
        }
        sites = predict_binding_sites(molecule)
        scored_sites = score_sites(sites)
        return {
            "sites": scored_sites,
            "count": len(scored_sites)
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Binding prediction failed: {str(e)}")


@router.post("/analysis")
async def analyze_route(req: MoleculeRequest):
    """Comprehensive molecular analysis"""
    try:
        molecule = {
            "atoms": req.atoms,
            "bonds": req.bonds
        }
        analysis = analyze_molecule(molecule)
        return analysis
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Analysis failed: {str(e)}")


@router.post("/analysis/similarity")
async def similarity_route(req: SimilarityRequest):
    """Compute similarity between two molecules"""
    try:
        mol1 = {
            "atoms": req.molecule1.atoms,
            "bonds": req.molecule1.bonds
        }
        mol2 = {
            "atoms": req.molecule2.atoms,
            "bonds": req.molecule2.bonds
        }
        similarity = compute_similarity(mol1, mol2)
        return {
            "similarity": round(similarity, 4),
            "similarity_percent": round(similarity * 100, 2)
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Similarity computation failed: {str(e)}")


@router.post("/analysis/alerts")
async def alerts_route(req: MoleculeRequest):
    """Get structural alerts for molecule"""
    try:
        molecule = {
            "atoms": req.atoms,
            "bonds": req.bonds
        }
        alerts = substructure_alerts(molecule)
        return {
            "alerts": alerts,
            "count": len(alerts)
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Alert detection failed: {str(e)}")


@router.post("/knowledge")
async def knowledge_route(req: MoleculeRequest):
    """Apply knowledge-based rules to molecule"""
    try:
        molecule = {
            "atoms": req.atoms,
            "bonds": req.bonds
        }
        rules = query_knowledge_rules(molecule)
        ml_predictions = apply_ml_rules(molecule)
        return {
            "rules": rules,
            "ml_predictions": ml_predictions,
            "rule_count": len(rules)
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Knowledge rules failed: {str(e)}")


@router.post("/explain")
async def explain_route(req: ExplainRequest):
    """Generate AI explanation for prediction/rule/alert"""
    try:
        molecule = {
            "atoms": req.molecule.atoms,
            "bonds": req.molecule.bonds
        }
        
        explanation = generate_explanation(molecule, req.site_id, req.alert_id, req.rule_id)
        return {
            "explanation": explanation
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Explanation generation failed: {str(e)}")


def generate_explanation(
    molecule: Dict[str, Any],
    site_id: Optional[int] = None,
    alert_id: Optional[int] = None,
    rule_id: Optional[str] = None
) -> str:
    """Generate textual explanation (simplified - can be enhanced with LLM)"""
    
    if site_id is not None:
        sites = predict_binding_sites(molecule)
        if 0 <= site_id < len(sites):
            site = sites[site_id]
            return (
                f"This binding site (type: {site.get('type', 'unknown')}) has a predicted "
                f"binding affinity score of {site.get('score', 0):.2f}. "
                f"The site is formed by {len(site.get('atom_indices', []))} atoms and is "
                f"characterized by {site.get('group_type', 'unknown')} functional groups. "
                f"Higher scores indicate stronger potential for ligand binding."
            )
        else:
            return "Binding site not found."
    
    if alert_id is not None:
        alerts = substructure_alerts(molecule)
        if 0 <= alert_id < len(alerts):
            alert = alerts[alert_id]
            return (
                f"Structural alert detected: {alert.get('type', 'unknown')} - "
                f"{alert.get('description', 'No description')}. "
                f"Severity: {alert.get('severity', 'unknown')}. "
                f"This alert affects {len(alert.get('affected_atoms', []))} atoms in the molecule."
            )
        else:
            return "Alert not found."
    
    if rule_id:
        rules = query_knowledge_rules(molecule)
        matching_rule = next((r for r in rules if r.get('rule') == rule_id), None)
        if matching_rule:
            return (
                f"Knowledge rule: {matching_rule.get('rule', 'unknown')} - "
                f"{matching_rule.get('description', 'No description')}. "
                f"Confidence: {matching_rule.get('confidence', 0):.2f}. "
                f"Category: {matching_rule.get('category', 'unknown')}."
            )
        else:
            return f"Rule '{rule_id}' not found."
    
    return "Please specify a site_id, alert_id, or rule_id to generate an explanation."

