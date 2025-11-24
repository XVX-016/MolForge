"""
Molecular Feature Extraction
Computes fingerprints and descriptors for ML models
"""
from typing import Dict, List, Any
from ..kab.utils import get_functional_groups, calculate_surface_area


def extract_features(molecule: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extract molecular features for ML prediction
    
    Args:
        molecule: Dict with 'atoms' and 'bonds' keys
        
    Returns:
        {
            "fingerprint": [0, 1, 0, 1, ...],  # Binary fingerprint
            "descriptors": {
                "molecular_weight": 180.16,
                "logp": 2.3,
                "hbd": 2,  # H-bond donors
                "hba": 4,  # H-bond acceptors
                "rotatable_bonds": 3,
                "aromatic_rings": 1
            }
        }
    """
    atoms = molecule.get('atoms', [])
    bonds = molecule.get('bonds', [])
    
    # Atomic weights
    atomic_weights = {
        'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999,
        'F': 18.998, 'P': 30.974, 'S': 32.065, 'Cl': 35.453,
        'Br': 79.904, 'I': 126.904
    }
    
    # Calculate molecular weight
    molecular_weight = sum(
        atomic_weights.get(atom.get('element', 'C'), 12.0)
        for atom in atoms
    )
    
    # Count H-bond donors (O-H, N-H)
    hbd = 0
    hba = 0
    
    # Build adjacency list
    adj = {i: [] for i in range(len(atoms))}
    for bond in bonds:
        a1 = bond.get('a1')
        a2 = bond.get('a2')
        if a1 is not None and a2 is not None:
            adj[a1].append(a2)
            adj[a2].append(a1)
    
    for i, atom in enumerate(atoms):
        element = atom.get('element', '')
        neighbors = adj.get(i, [])
        
        # H-bond donors
        if element in ['O', 'N']:
            for neighbor_idx in neighbors:
                neighbor = atoms[neighbor_idx] if neighbor_idx < len(atoms) else None
                if neighbor and neighbor.get('element') == 'H':
                    hbd += 1
                    break
        
        # H-bond acceptors (O, N with lone pairs)
        if element in ['O', 'N']:
            max_bonds = 2 if element == 'O' else 3
            if len(neighbors) < max_bonds:
                hba += 1
    
    # Count rotatable bonds (single bonds not in rings)
    rotatable_bonds = 0
    for bond in bonds:
        if bond.get('order', 1) == 1:
            # Simplified: assume not in ring (would need ring detection)
            rotatable_bonds += 1
    
    # Count aromatic rings (simplified)
    functional_groups = get_functional_groups(molecule)
    aromatic_rings = sum(1 for g in functional_groups if g.get('type') == 'aromatic')
    
    # Estimate LogP (simplified: based on MW and polarity)
    # LogP ≈ 0.5 * (non-polar atoms) - 0.5 * (polar atoms)
    non_polar = sum(1 for a in atoms if a.get('element') in ['C', 'H'])
    polar = sum(1 for a in atoms if a.get('element') in ['O', 'N', 'F', 'Cl'])
    logp = 0.5 * non_polar - 0.5 * polar
    
    # Simple binary fingerprint (MACCS-like, simplified)
    fingerprint = generate_simple_fingerprint(molecule)
    
    return {
        "fingerprint": fingerprint,
        "descriptors": {
            "molecular_weight": round(molecular_weight, 2),
            "logp": round(logp, 2),
            "hbd": hbd,
            "hba": hba,
            "rotatable_bonds": rotatable_bonds,
            "aromatic_rings": aromatic_rings,
            "num_atoms": len(atoms),
            "num_bonds": len(bonds)
        }
    }


def generate_simple_fingerprint(molecule: Dict[str, Any]) -> List[int]:
    """
    Generate simple binary fingerprint (MACCS-like)
    
    Returns 166-bit fingerprint (simplified to 32 bits for demo)
    """
    atoms = molecule.get('atoms', [])
    bonds = molecule.get('bonds', [])
    functional_groups = get_functional_groups(molecule)
    
    # 32-bit fingerprint
    fp = [0] * 32
    
    # Bit 0: Has carbon
    if any(a.get('element') == 'C' for a in atoms):
        fp[0] = 1
    
    # Bit 1: Has nitrogen
    if any(a.get('element') == 'N' for a in atoms):
        fp[1] = 1
    
    # Bit 2: Has oxygen
    if any(a.get('element') == 'O' for a in atoms):
        fp[2] = 1
    
    # Bit 3: Has aromatic ring
    if any(g.get('type') == 'aromatic' for g in functional_groups):
        fp[3] = 1
    
    # Bit 4: Has carbonyl
    if any(g.get('type') == 'carbonyl' for g in functional_groups):
        fp[4] = 1
    
    # Bit 5: Has hydroxyl
    if any(g.get('type') == 'hydroxyl' for g in functional_groups):
        fp[5] = 1
    
    # Bit 6: Has amino
    if any(g.get('type') == 'amino' for g in functional_groups):
        fp[6] = 1
    
    # Bit 7: Has double bonds
    if any(b.get('order', 1) >= 2 for b in bonds):
        fp[7] = 1
    
    # Remaining bits: element counts and other features
    element_counts = {}
    for atom in atoms:
        element = atom.get('element', 'C')
        element_counts[element] = element_counts.get(element, 0) + 1
    
    # Bits 8-15: Element presence flags
    elements = ['H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl']
    for i, elem in enumerate(elements):
        if i + 8 < 32 and elem in element_counts:
            fp[i + 8] = 1
    
    # Bits 16-23: Size features
    if len(atoms) > 10:
        fp[16] = 1
    if len(atoms) > 20:
        fp[17] = 1
    if len(bonds) > 10:
        fp[18] = 1
    
    return fp

