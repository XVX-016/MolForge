"""
ML Property Prediction Engine
Uses extracted features to predict molecular properties
"""
from typing import Dict, List, Any
from .features import extract_features


def predict_properties(molecule: Dict[str, Any]) -> Dict[str, Any]:
    """
    Predict molecular properties using ML models (or heuristics)
    
    Args:
        molecule: Dict with 'atoms' and 'bonds' keys
        
    Returns:
        {
            "properties": [
                {"property": "solubility", "value": 0.78, "unit": "logS"},
                {"property": "toxicity", "value": 0.23, "unit": "probability"},
                {"property": "bioavailability", "value": 0.65, "unit": "fraction"}
            ],
            "model": "heuristic"
        }
    """
    features = extract_features(molecule)
    descriptors = features.get('descriptors', {})
    
    # Heuristic-based predictions (can be replaced with actual ML models)
    predictions = []
    
    # Solubility (logS) - simplified model
    mw = descriptors.get('molecular_weight', 200)
    logp = descriptors.get('logp', 2.0)
    hbd = descriptors.get('hbd', 0)
    hba = descriptors.get('hba', 0)
    
    # Simplified: solubility increases with H-bonding, decreases with MW and LogP
    solubility = 0.5 - (logp * 0.1) - (mw / 1000) + (hbd + hba) * 0.05
    solubility = max(-6.0, min(1.0, solubility))  # Clamp to reasonable range
    predictions.append({
        "property": "solubility",
        "value": round(solubility, 2),
        "unit": "logS"
    })
    
    # Toxicity risk (simplified)
    # Higher risk for heavy metals, reactive groups
    toxicity = 0.2  # Base risk
    atoms = molecule.get('atoms', [])
    for atom in atoms:
        element = atom.get('element', '')
        if element in ['Pb', 'Hg', 'Cd', 'As']:
            toxicity += 0.3
        elif element in ['Br', 'I']:
            toxicity += 0.1
    
    toxicity = min(1.0, toxicity)
    predictions.append({
        "property": "toxicity",
        "value": round(toxicity, 2),
        "unit": "probability"
    })
    
    # Bioavailability (fraction)
    # Based on Lipinski's Rule of Five
    bioavailability = 0.7  # Base
    if mw > 500:
        bioavailability -= 0.2
    if hbd > 5:
        bioavailability -= 0.1
    if hba > 10:
        bioavailability -= 0.1
    if logp > 5:
        bioavailability -= 0.1
    
    bioavailability = max(0.0, min(1.0, bioavailability))
    predictions.append({
        "property": "bioavailability",
        "value": round(bioavailability, 2),
        "unit": "fraction"
    })
    
    # Drug-likeness score
    drug_likeness = 0.5
    if 200 <= mw <= 500:
        drug_likeness += 0.1
    if hbd <= 5:
        drug_likeness += 0.1
    if hba <= 10:
        drug_likeness += 0.1
    if -2 <= logp <= 5:
        drug_likeness += 0.2
    
    drug_likeness = min(1.0, drug_likeness)
    predictions.append({
        "property": "drug_likeness",
        "value": round(drug_likeness, 2),
        "unit": "score"
    })
    
    return {
        "properties": predictions,
        "model": "heuristic",
        "features": descriptors
    }

