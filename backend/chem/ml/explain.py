"""
AI Explainability Engine
Provides per-atom/per-bond explanations for ML predictions and retrosynthesis
"""
from typing import Dict, List, Any, Optional
from .features import extract_features
from .predict import predict_properties
from ..kab.utils import get_atom_neighbors, get_functional_groups


def explain_prediction(
    molecule: Dict[str, Any],
    property_name: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate explanation for ML property prediction
    
    Args:
        molecule: Dict with 'atoms' and 'bonds' keys
        property_name: Optional specific property to explain (e.g., 'solubility')
        
    Returns:
        {
            "property": "solubility",
            "value": 0.78,
            "explanation": "Hydrophilic groups on atoms 2,3 increase solubility",
            "atom_contributions": [
                {"atom": 2, "contribution": 0.12, "reason": "Hydroxyl group"},
                {"atom": 3, "contribution": 0.15, "reason": "Amino group"}
            ],
            "bond_contributions": [...]
        }
    """
    atoms = molecule.get('atoms', [])
    bonds = molecule.get('bonds', [])
    
    if not atoms:
        return {
            "property": property_name or "unknown",
            "value": 0.0,
            "explanation": "No molecule data available",
            "atom_contributions": [],
            "bond_contributions": []
        }
    
    # Get predictions
    predictions = predict_properties(molecule)
    properties = predictions.get('properties', [])
    
    # If specific property requested, find it
    target_prop = None
    if property_name:
        target_prop = next((p for p in properties if p.get('property') == property_name), None)
    else:
        # Default to first property
        target_prop = properties[0] if properties else None
    
    if not target_prop:
        return {
            "property": property_name or "unknown",
            "value": 0.0,
            "explanation": "Property not found",
            "atom_contributions": [],
            "bond_contributions": []
        }
    
    prop_name = target_prop.get('property', 'unknown')
    prop_value = target_prop.get('value', 0.0)
    
    # Calculate atom contributions using feature importance
    atom_contributions = calculate_atom_contributions(molecule, prop_name, prop_value)
    
    # Generate textual explanation
    explanation = generate_textual_explanation(molecule, prop_name, prop_value, atom_contributions)
    
    return {
        "property": prop_name,
        "value": prop_value,
        "unit": target_prop.get('unit', ''),
        "explanation": explanation,
        "atom_contributions": atom_contributions,
        "bond_contributions": []  # Placeholder for bond contributions
    }


def calculate_atom_contributions(
    molecule: Dict[str, Any],
    property_name: str,
    property_value: float
) -> List[Dict[str, Any]]:
    """
    Calculate per-atom contributions to a property prediction
    
    Uses simplified feature importance based on functional groups and local environment
    """
    atoms = molecule.get('atoms', [])
    bonds = molecule.get('bonds', [])
    functional_groups = get_functional_groups(molecule)
    
    if not atoms:
        return []
    
    # Build adjacency list
    adj = {i: [] for i in range(len(atoms))}
    for bond in bonds:
        a1 = bond.get('a1')
        a2 = bond.get('a2')
        if a1 is not None and a2 is not None:
            adj[a1].append(a2)
            adj[a2].append(a1)
    
    contributions = []
    
    for i, atom in enumerate(atoms):
        element = atom.get('element', 'C')
        neighbors = adj.get(i, [])
        
        contribution = 0.0
        reasons = []
        
        # Property-specific contribution rules
        if property_name == 'solubility':
            # Hydrophilic groups increase solubility
            if element in ['O', 'N']:
                # Check if part of polar functional group
                for group in functional_groups:
                    if i in group.get('atoms', []):
                        group_type = group.get('type', '')
                        if group_type in ['hydroxyl', 'amino', 'carbonyl']:
                            contribution += 0.15
                            reasons.append(f"{group_type.capitalize()} group")
                            break
                
                # Lone pairs (H-bond acceptors)
                max_bonds = 2 if element == 'O' else 3
                if len(neighbors) < max_bonds:
                    contribution += 0.08
                    reasons.append("H-bond acceptor")
            
            # Hydrophobic groups decrease solubility
            if element == 'C' and len(neighbors) >= 3:
                # Check if in hydrophobic region
                carbon_neighbors = sum(1 for n in neighbors if atoms[n].get('element') == 'C')
                if carbon_neighbors >= 2:
                    contribution -= 0.05
                    reasons.append("Hydrophobic carbon chain")
        
        elif property_name == 'toxicity':
            # Heavy metals and reactive groups increase toxicity
            if element in ['Pb', 'Hg', 'Cd', 'As']:
                contribution += 0.4
                reasons.append(f"Heavy metal ({element})")
            elif element in ['Br', 'I']:
                contribution += 0.1
                reasons.append(f"Halogen ({element})")
            
            # Check for reactive groups
            for group in functional_groups:
                if i in group.get('atoms', []):
                    group_type = group.get('type', '')
                    if group_type in ['peroxide', 'epoxide']:
                        contribution += 0.2
                        reasons.append("Reactive group")
                        break
        
        elif property_name == 'bioavailability':
            # Polar groups can help or hinder depending on context
            if element in ['O', 'N']:
                # H-bond donors/acceptors can help
                for neighbor_idx in neighbors:
                    neighbor = atoms[neighbor_idx] if neighbor_idx < len(atoms) else None
                    if neighbor and neighbor.get('element') == 'H':
                        contribution += 0.1
                        reasons.append("H-bond donor")
                        break
        
        elif property_name == 'drug_likeness':
            # Lipinski's Rule of Five factors
            if element in ['O', 'N']:
                # Too many H-bond donors/acceptors can hurt
                if len(neighbors) < (2 if element == 'O' else 3):
                    contribution -= 0.05
                    reasons.append("H-bond acceptor")
        
        # Normalize contribution
        contribution = max(-0.5, min(0.5, contribution))
        
        if abs(contribution) > 0.01:  # Only include significant contributions
            contributions.append({
                "atom": i,
                "contribution": round(contribution, 3),
                "reason": "; ".join(reasons) if reasons else "Local environment",
                "element": element
            })
    
    # Sort by absolute contribution (most important first)
    contributions.sort(key=lambda x: abs(x['contribution']), reverse=True)
    
    return contributions


def generate_textual_explanation(
    molecule: Dict[str, Any],
    property_name: str,
    property_value: float,
    atom_contributions: List[Dict[str, Any]]
) -> str:
    """
    Generate human-readable explanation for property prediction
    """
    atoms = molecule.get('atoms', [])
    
    # Get top contributing atoms
    top_contributors = atom_contributions[:5]  # Top 5 contributors
    
    if not top_contributors:
        return f"The predicted {property_name} value of {property_value:.2f} is based on overall molecular structure."
    
    positive_contributors = [c for c in top_contributors if c['contribution'] > 0]
    negative_contributors = [c for c in top_contributors if c['contribution'] < 0]
    
    explanation_parts = []
    
    if positive_contributors:
        atom_indices = [str(c['atom']) for c in positive_contributors]
        reasons = [c['reason'] for c in positive_contributors if c.get('reason')]
        explanation_parts.append(
            f"Atoms {', '.join(atom_indices)} contribute positively to {property_name} "
            f"({', '.join(set(reasons[:3])) if reasons else 'structural features'})."
        )
    
    if negative_contributors:
        atom_indices = [str(c['atom']) for c in negative_contributors]
        reasons = [c['reason'] for c in negative_contributors if c.get('reason')]
        explanation_parts.append(
            f"Atoms {', '.join(atom_indices)} contribute negatively to {property_name} "
            f"({', '.join(set(reasons[:3])) if reasons else 'structural features'})."
        )
    
    if not explanation_parts:
        explanation_parts.append(
            f"The predicted {property_name} value of {property_value:.2f} is based on overall molecular structure."
        )
    
    return " ".join(explanation_parts)


def explain_retrosynthesis_step(
    pathway_step: Dict[str, Any],
    reaction: Optional[Dict[str, Any]] = None
) -> str:
    """
    Generate explanation for a retrosynthesis step
    
    Args:
        pathway_step: Step data from retrosynthesis planner
        reaction: Optional reaction data
        
    Returns:
        Human-readable explanation string
    """
    step_num = pathway_step.get('step', 0)
    precursors = pathway_step.get('precursors', [])
    reaction_name = reaction.get('name', 'Unknown reaction') if reaction else 'Unknown reaction'
    
    if step_num == 0:
        return f"Starting material: {len(precursors)} precursor(s) identified."
    
    explanation = f"Step {step_num}: Apply {reaction_name} to generate "
    
    if precursors:
        explanation += f"{len(precursors)} precursor molecule(s). "
        if len(precursors) == 1:
            precursor = precursors[0]
            atoms = precursor.get('atoms', [])
            explanation += f"The precursor contains {len(atoms)} atoms."
        else:
            explanation += "These precursors can be further deconstructed."
    else:
        explanation += "No valid precursors found."
    
    return explanation


def explain_reaction_mechanism(
    mechanism_step: Dict[str, Any],
    step_index: int
) -> str:
    """
    Generate explanation for a reaction mechanism step
    
    Args:
        mechanism_step: Step data from mechanism prediction
        step_index: Index of step in mechanism
        
    Returns:
        Human-readable explanation string
    """
    description = mechanism_step.get('description', f'Step {step_index}')
    energy = mechanism_step.get('energy', 0.0)
    
    explanation = f"{description}. "
    
    if step_index == 0:
        explanation += "This is the initial reactant state."
    else:
        if energy > 0:
            explanation += f"The intermediate has an energy of {energy:.2f} kcal/mol above the reactant, indicating an energy barrier."
        else:
            explanation += f"The intermediate has an energy of {energy:.2f} kcal/mol, lower than the reactant, indicating a favorable step."
    
    return explanation

