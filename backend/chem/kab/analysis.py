"""
Molecular Analysis & Similarity Engine
Computes similarity, substructure analysis, and molecular descriptors
"""
from typing import Dict, List, Any, Optional, Tuple
from .utils import get_functional_groups, calculate_surface_area


def compute_similarity(mol1: Dict[str, Any], mol2: Dict[str, Any]) -> float:
    """
    Compute Tanimoto similarity between two molecules
    
    Uses functional group fingerprints for similarity calculation
    
    Args:
        mol1: First molecule dict
        mol2: Second molecule dict
        
    Returns:
        Similarity score (0-1)
    """
    # Get functional groups for both molecules
    groups1 = get_functional_groups(mol1)
    groups2 = get_functional_groups(mol2)
    
    # Create fingerprints based on functional group types
    fp1 = set(g.get('type', '') for g in groups1)
    fp2 = set(g.get('type', '') for g in groups2)
    
    # Also include element composition
    atoms1 = mol1.get('atoms', [])
    atoms2 = mol2.get('atoms', [])
    elements1 = set(atom.get('element', '') for atom in atoms1)
    elements2 = set(atom.get('element', '') for atom in atoms2)
    
    # Combine fingerprints
    fp1_full = fp1.union(elements1)
    fp2_full = fp2.union(elements2)
    
    # Tanimoto coefficient: intersection / union
    intersection = len(fp1_full.intersection(fp2_full))
    union = len(fp1_full.union(fp2_full))
    
    if union == 0:
        return 0.0
    
    similarity = intersection / union
    return similarity


def substructure_alerts(molecule: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Detect structural alerts (problematic substructures)
    
    Args:
        molecule: Molecule dict
        
    Returns:
        List of alerts with:
        - type: Alert type (e.g., 'toxic', 'reactive', 'unstable')
        - description: Human-readable description
        - affected_atoms: List of atom indices
        - severity: 'low', 'medium', 'high'
    """
    atoms = molecule.get('atoms', [])
    bonds = molecule.get('bonds', [])
    alerts = []
    
    # Build adjacency list
    adj = {i: [] for i in range(len(atoms))}
    for bond in bonds:
        a1 = bond.get('a1')
        a2 = bond.get('a2')
        if a1 is not None and a2 is not None:
            adj[a1].append(a2)
            adj[a2].append(a1)
    
    # Check for problematic patterns
    for i, atom in enumerate(atoms):
        element = atom.get('element', '')
        neighbors = adj.get(i, [])
        
        # Highly reactive: epoxide (3-membered ring with O)
        if element == 'O' and len(neighbors) == 2:
            # Check if forms small ring
            n1_idx = neighbors[0]
            n2_idx = neighbors[1]
            n1_neighbors = adj.get(n1_idx, [])
            n2_neighbors = adj.get(n2_idx, [])
            
            # If neighbors are connected, might be epoxide
            if n2_idx in n1_neighbors:
                alerts.append({
                    'type': 'reactive',
                    'description': 'Potential epoxide ring detected (highly reactive)',
                    'affected_atoms': [i, n1_idx, n2_idx],
                    'severity': 'high'
                })
        
        # Unstable: peroxide (O-O bond)
        if element == 'O':
            for neighbor_idx in neighbors:
                neighbor = atoms[neighbor_idx] if neighbor_idx < len(atoms) else None
                if neighbor and neighbor.get('element') == 'O':
                    alerts.append({
                        'type': 'unstable',
                        'description': 'Peroxide bond detected (O-O, potentially explosive)',
                        'affected_atoms': [i, neighbor_idx],
                        'severity': 'high'
                    })
        
        # Toxic: heavy metals (simplified check)
        toxic_elements = ['Pb', 'Hg', 'Cd', 'As', 'Cr']
        if element in toxic_elements:
            alerts.append({
                'type': 'toxic',
                'description': f'Heavy metal {element} detected (potential toxicity)',
                'affected_atoms': [i],
                'severity': 'high'
            })
    
    # Check for strained rings (small rings are often unstable)
    ring_sizes = detect_rings(molecule)
    for ring in ring_sizes:
        if ring['size'] == 3:
            alerts.append({
                'type': 'unstable',
                'description': '3-membered ring detected (highly strained, unstable)',
                'affected_atoms': ring['atoms'],
                'severity': 'medium'
            })
        elif ring['size'] == 4:
            alerts.append({
                'type': 'unstable',
                'description': '4-membered ring detected (strained, may be unstable)',
                'affected_atoms': ring['atoms'],
                'severity': 'low'
            })
    
    return alerts


def detect_rings(molecule: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Detect rings in molecule (simplified DFS-based detection)"""
    atoms = molecule.get('atoms', [])
    bonds = molecule.get('bonds', [])
    
    if len(atoms) < 3:
        return []
    
    # Build adjacency list
    adj = {i: [] for i in range(len(atoms))}
    for bond in bonds:
        a1 = bond.get('a1')
        a2 = bond.get('a2')
        if a1 is not None and a2 is not None:
            adj[a1].append(a2)
            adj[a2].append(a1)
    
    rings = []
    visited = set()
    
    def dfs_ring(start: int, current: int, path: List[int], min_ring_size: int = 3, max_ring_size: int = 8):
        """DFS to find rings"""
        if len(path) > max_ring_size:
            return
        
        for neighbor in adj.get(current, []):
            if neighbor == start and len(path) >= min_ring_size:
                # Found a ring
                ring_atoms = path + [start]
                rings.append({
                    'atoms': ring_atoms,
                    'size': len(ring_atoms)
                })
                return
            elif neighbor not in path and neighbor not in visited:
                dfs_ring(start, neighbor, path + [neighbor], min_ring_size, max_ring_size)
    
    # Try to find rings starting from each atom
    for start in range(len(atoms)):
        if start not in visited:
            dfs_ring(start, start, [start])
            visited.add(start)
    
    # Remove duplicates (same ring, different starting points)
    unique_rings = []
    seen_rings = set()
    for ring in rings:
        ring_set = tuple(sorted(ring['atoms']))
        if ring_set not in seen_rings:
            seen_rings.add(ring_set)
            unique_rings.append(ring)
    
    return unique_rings


def analyze_molecule(molecule: Dict[str, Any]) -> Dict[str, Any]:
    """
    Comprehensive molecular analysis
    
    Args:
        molecule: Molecule dict
        
    Returns:
        Analysis results including:
        - functional_groups: List of detected groups
        - molecular_weight: Estimated MW
        - surface_area: Estimated surface area
        - complexity: Molecular complexity score
        - alerts: Structural alerts
    """
    atoms = molecule.get('atoms', [])
    bonds = molecule.get('bonds', [])
    
    # Functional groups
    functional_groups = get_functional_groups(molecule)
    
    # Molecular weight (simplified)
    atomic_weights = {
        'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999,
        'F': 18.998, 'P': 30.974, 'S': 32.065, 'Cl': 35.453,
        'Br': 79.904, 'I': 126.904
    }
    
    molecular_weight = sum(
        atomic_weights.get(atom.get('element', 'C'), 12.0)
        for atom in atoms
    )
    
    # Surface area
    surface_area = calculate_surface_area(molecule)
    
    # Complexity score (based on number of atoms, bonds, rings)
    rings = detect_rings(molecule)
    complexity = (
        len(atoms) * 0.1 +
        len(bonds) * 0.05 +
        len(rings) * 0.2 +
        len(functional_groups) * 0.15
    )
    
    # Structural alerts
    alerts = substructure_alerts(molecule)
    
    return {
        'functional_groups': functional_groups,
        'molecular_weight': round(molecular_weight, 2),
        'surface_area': round(surface_area, 2),
        'complexity': round(complexity, 2),
        'alerts': alerts,
        'num_atoms': len(atoms),
        'num_bonds': len(bonds),
        'num_rings': len(rings)
    }

