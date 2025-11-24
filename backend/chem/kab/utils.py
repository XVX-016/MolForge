"""
Utility functions for KAB engine
"""
from typing import Dict, List, Any, Tuple
import math


def calculate_distance(pos1: List[float], pos2: List[float]) -> float:
    """Calculate Euclidean distance between two 3D points"""
    return math.sqrt(
        (pos1[0] - pos2[0]) ** 2 +
        (pos1[1] - pos2[1]) ** 2 +
        (pos1[2] - pos2[2]) ** 2
    )


def get_atom_neighbors(
    molecule: Dict[str, Any],
    atom_index: int,
    max_distance: float = 3.0
) -> List[int]:
    """Get neighboring atoms within max_distance"""
    atoms = molecule.get('atoms', [])
    bonds = molecule.get('bonds', [])
    
    if atom_index >= len(atoms):
        return []
    
    atom = atoms[atom_index]
    neighbors = []
    
    # First, check direct bonds
    for bond in bonds:
        if bond.get('a1') == atom_index:
            neighbors.append(bond.get('a2'))
        elif bond.get('a2') == atom_index:
            neighbors.append(bond.get('a1'))
    
    # Then check spatial neighbors
    atom_pos = atom.get('position', [0, 0, 0])
    for i, other_atom in enumerate(atoms):
        if i == atom_index:
            continue
        if i in neighbors:
            continue
        
        other_pos = other_atom.get('position', [0, 0, 0])
        distance = calculate_distance(atom_pos, other_pos)
        if distance <= max_distance:
            neighbors.append(i)
    
    return neighbors


def get_functional_groups(molecule: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Detect functional groups in molecule"""
    atoms = molecule.get('atoms', [])
    bonds = molecule.get('bonds', [])
    groups = []
    
    # Build adjacency list
    adj = {i: [] for i in range(len(atoms))}
    for bond in bonds:
        a1 = bond.get('a1')
        a2 = bond.get('a2')
        if a1 is not None and a2 is not None:
            adj[a1].append(a2)
            adj[a2].append(a1)
    
    # Detect common functional groups
    for i, atom in enumerate(atoms):
        element = atom.get('element', '')
        
        # Hydroxyl group (OH)
        if element == 'O':
            neighbors = adj.get(i, [])
            if len(neighbors) == 1:  # Terminal oxygen
                neighbor = atoms[neighbors[0]] if neighbors else None
                if neighbor and neighbor.get('element') == 'H':
                    groups.append({
                        'type': 'hydroxyl',
                        'atoms': [i, neighbors[0]],
                        'center': i
                    })
        
        # Carbonyl group (C=O)
        if element == 'C':
            neighbors = adj.get(i, [])
            for neighbor_idx in neighbors:
                neighbor = atoms[neighbor_idx] if neighbor_idx < len(atoms) else None
                if neighbor and neighbor.get('element') == 'O':
                    # Check if double bond
                    bond_order = 1
                    for bond in bonds:
                        if (bond.get('a1') == i and bond.get('a2') == neighbor_idx) or \
                           (bond.get('a1') == neighbor_idx and bond.get('a2') == i):
                            bond_order = bond.get('order', 1)
                            break
                    if bond_order == 2:
                        groups.append({
                            'type': 'carbonyl',
                            'atoms': [i, neighbor_idx],
                            'center': i
                        })
        
        # Amino group (NH2)
        if element == 'N':
            neighbors = adj.get(i, [])
            h_count = sum(1 for n in neighbors if atoms[n].get('element') == 'H')
            if h_count >= 2:
                h_atoms = [n for n in neighbors if atoms[n].get('element') == 'H']
                groups.append({
                    'type': 'amino',
                    'atoms': [i] + h_atoms[:2],
                    'center': i
                })
    
    return groups


def calculate_surface_area(molecule: Dict[str, Any]) -> float:
    """Estimate molecular surface area (simplified)"""
    atoms = molecule.get('atoms', [])
    if not atoms:
        return 0.0
    
    # Van der Waals radii (approximate, in Angstroms)
    vdw_radii = {
        'H': 1.2, 'C': 1.7, 'N': 1.55, 'O': 1.52,
        'F': 1.47, 'P': 1.8, 'S': 1.8, 'Cl': 1.75,
        'Br': 1.85, 'I': 1.98
    }
    
    total_area = 0.0
    for atom in atoms:
        element = atom.get('element', 'C')
        radius = vdw_radii.get(element, 1.7)
        # Surface area of sphere: 4πr²
        area = 4 * math.pi * (radius ** 2)
        total_area += area
    
    # Rough correction for overlapping spheres
    return total_area * 0.7  # Approximate reduction factor

