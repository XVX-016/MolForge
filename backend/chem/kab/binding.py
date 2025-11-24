"""
Binding Site Prediction Engine
Predicts potential ligand binding sites based on molecular structure
"""
from typing import Dict, List, Any
from .utils import get_atom_neighbors, get_functional_groups, calculate_surface_area


def predict_binding_sites(molecule: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Predict potential binding sites in molecule
    
    Args:
        molecule: Dict with 'atoms' and 'bonds' keys
        
    Returns:
        List of binding sites with:
        - atom_indices: List of atom indices forming the site
        - score: Binding affinity score (0-1)
        - type: Type of binding site (e.g., 'hydrophobic', 'polar', 'aromatic')
        - center: Center position of the site
    """
    atoms = molecule.get('atoms', [])
    if not atoms:
        return []
    
    sites = []
    
    # Get functional groups (potential binding sites)
    functional_groups = get_functional_groups(molecule)
    
    # Score each functional group as a potential binding site
    for group in functional_groups:
        group_type = group.get('type', '')
        group_atoms = group.get('atoms', [])
        
        if not group_atoms:
            continue
        
        # Calculate binding score based on group type and environment
        score = calculate_binding_score(molecule, group)
        
        # Get center position
        center_atom_idx = group.get('center', group_atoms[0])
        center_atom = atoms[center_atom_idx] if center_atom_idx < len(atoms) else None
        center = center_atom.get('position', [0, 0, 0]) if center_atom else [0, 0, 0]
        
        # Determine site type
        site_type = map_group_to_site_type(group_type)
        
        sites.append({
            'atom_indices': group_atoms,
            'score': score,
            'type': site_type,
            'center': center,
            'group_type': group_type
        })
    
    # Also identify hydrophobic pockets
    hydrophobic_sites = identify_hydrophobic_pockets(molecule)
    sites.extend(hydrophobic_sites)
    
    # Sort by score (highest first)
    sites.sort(key=lambda x: x['score'], reverse=True)
    
    return sites


def calculate_binding_score(molecule: Dict[str, Any], group: Dict[str, Any]) -> float:
    """Calculate binding affinity score for a functional group"""
    group_type = group.get('type', '')
    atoms = molecule.get('atoms', [])
    group_atoms = group.get('atoms', [])
    
    base_scores = {
        'hydroxyl': 0.7,      # Polar, can form H-bonds
        'carbonyl': 0.8,      # Strong polar interaction
        'amino': 0.75,        # Can form H-bonds, charged
        'carboxyl': 0.85,    # Strong polar, can be charged
        'aromatic': 0.6,     # π-π stacking
    }
    
    base_score = base_scores.get(group_type, 0.5)
    
    # Adjust based on accessibility (surface exposure)
    # Simplified: check if group has many neighbors (buried) or few (exposed)
    center_idx = group.get('center', group_atoms[0] if group_atoms else 0)
    neighbors = get_atom_neighbors(molecule, center_idx, max_distance=4.0)
    
    # More exposed = higher score
    exposure_factor = 1.0 - (len(neighbors) / 10.0)  # Normalize
    exposure_factor = max(0.3, min(1.0, exposure_factor))
    
    # Adjust based on electronegativity (polar groups score higher)
    if group_type in ['hydroxyl', 'carbonyl', 'amino']:
        polarity_bonus = 0.1
    else:
        polarity_bonus = 0.0
    
    final_score = base_score * exposure_factor + polarity_bonus
    return min(1.0, max(0.0, final_score))


def map_group_to_site_type(group_type: str) -> str:
    """Map functional group type to binding site type"""
    mapping = {
        'hydroxyl': 'polar',
        'carbonyl': 'polar',
        'amino': 'polar',
        'carboxyl': 'polar',
        'aromatic': 'aromatic',
    }
    return mapping.get(group_type, 'hydrophobic')


def identify_hydrophobic_pockets(molecule: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Identify hydrophobic regions (carbon-rich areas)"""
    atoms = molecule.get('atoms', [])
    if not atoms:
        return []
    
    sites = []
    
    # Find clusters of carbon atoms
    carbon_atoms = [i for i, atom in enumerate(atoms) if atom.get('element') == 'C']
    
    if len(carbon_atoms) < 3:
        return sites
    
    # Group nearby carbons
    clusters = []
    for carbon_idx in carbon_atoms:
        # Check if this carbon is near an existing cluster
        added = False
        for cluster in clusters:
            for cluster_idx in cluster:
                neighbors = get_atom_neighbors(molecule, carbon_idx, max_distance=3.0)
                if cluster_idx in neighbors:
                    cluster.append(carbon_idx)
                    added = True
                    break
            if added:
                break
        
        if not added:
            clusters.append([carbon_idx])
    
    # Create binding sites for large carbon clusters
    for cluster in clusters:
        if len(cluster) >= 3:  # At least 3 carbons
            # Calculate center
            positions = [atoms[i].get('position', [0, 0, 0]) for i in cluster]
            center = [
                sum(p[0] for p in positions) / len(positions),
                sum(p[1] for p in positions) / len(positions),
                sum(p[2] for p in positions) / len(positions),
            ]
            
            # Score based on cluster size
            score = min(0.6, 0.3 + (len(cluster) / 10.0))
            
            sites.append({
                'atom_indices': cluster,
                'score': score,
                'type': 'hydrophobic',
                'center': center,
                'group_type': 'hydrophobic_pocket'
            })
    
    return sites


def score_sites(sites: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Normalize and score binding sites
    
    Args:
        sites: List of binding sites from predict_binding_sites
        
    Returns:
        Sites with normalized scores (0-1)
    """
    if not sites:
        return []
    
    # Get max score for normalization
    max_score = max(site.get('score', 0) for site in sites) if sites else 1.0
    
    if max_score == 0:
        max_score = 1.0
    
    # Normalize scores
    normalized = []
    for site in sites:
        normalized_site = site.copy()
        normalized_site['score'] = site.get('score', 0) / max_score
        normalized_site['normalized_score'] = normalized_site['score']
        normalized.append(normalized_site)
    
    return normalized

