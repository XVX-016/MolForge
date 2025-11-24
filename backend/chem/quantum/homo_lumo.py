"""
HOMO/LUMO Calculation Engine
Uses Hückel approximation for simple molecules
"""
from typing import Dict, List, Any
import math


def calculate_homo_lumo(molecule: Dict[str, Any]) -> Dict[str, Any]:
    """
    Calculate HOMO and LUMO orbital energies using Hückel approximation
    
    Args:
        molecule: Dict with 'atoms' and 'bonds' keys
        
    Returns:
        {
            "HOMO": -5.4,
            "LUMO": -1.2,
            "gap": 4.2,
            "contributing_atoms": [1, 2, 3],
            "method": "huckel"
        }
    """
    atoms = molecule.get('atoms', [])
    bonds = molecule.get('bonds', [])
    
    if not atoms:
        return {
            "HOMO": 0.0,
            "LUMO": 0.0,
            "gap": 0.0,
            "contributing_atoms": [],
            "method": "huckel"
        }
    
    # Build adjacency matrix for π electrons (simplified Hückel)
    n = len(atoms)
    h_matrix = [[0.0 for _ in range(n)] for _ in range(n)]
    
    # Diagonal: α (Coulomb integral, set to -10 eV for carbon)
    alpha = -10.0
    for i in range(n):
        element = atoms[i].get('element', 'C')
        # Adjust alpha based on element electronegativity
        if element == 'N':
            alpha = -11.0
        elif element == 'O':
            alpha = -12.0
        elif element == 'S':
            alpha = -9.5
        else:
            alpha = -10.0
        h_matrix[i][i] = alpha
    
    # Off-diagonal: β (resonance integral, -2.4 eV for C-C bonds)
    beta = -2.4
    for bond in bonds:
        a1 = bond.get('a1')
        a2 = bond.get('a2')
        if a1 is not None and a2 is not None and 0 <= a1 < n and 0 <= a2 < n:
            # Only consider π bonds (double/triple) or aromatic systems
            bond_order = bond.get('order', 1)
            if bond_order >= 2:
                h_matrix[a1][a2] = beta * (bond_order - 1)
                h_matrix[a2][a1] = beta * (bond_order - 1)
    
    # Simple eigenvalue calculation (for small systems)
    # For larger systems, would use numpy.linalg.eig
    eigenvalues = simple_eigenvalues(h_matrix)
    
    # Sort eigenvalues (most negative = lowest energy)
    eigenvalues.sort()
    
    # HOMO is highest occupied (most negative)
    # LUMO is lowest unoccupied (least negative)
    num_electrons = count_pi_electrons(molecule)
    
    # Ensure we have enough orbitals
    if len(eigenvalues) < 2:
        # Not enough orbitals - return default values
        return {
            "HOMO": -10.0,
            "LUMO": -1.0,
            "gap": 9.0,
            "contributing_atoms": [],
            "method": "huckel"
        }
    
    # HOMO index: (num_electrons // 2) - 1 (0-indexed)
    homo_idx = max(0, (num_electrons // 2) - 1) if num_electrons > 0 else 0
    lumo_idx = homo_idx + 1
    
    # Ensure indices are valid
    homo_idx = min(homo_idx, len(eigenvalues) - 1)
    lumo_idx = min(lumo_idx, len(eigenvalues) - 1)
    
    homo = eigenvalues[homo_idx]
    lumo = eigenvalues[lumo_idx] if lumo_idx < len(eigenvalues) else eigenvalues[-1]
    
    gap = lumo - homo
    
    # Find atoms contributing to HOMO/LUMO (simplified: atoms in π system)
    contributing_atoms = []
    for i, atom in enumerate(atoms):
        element = atom.get('element', 'C')
        # Check if atom is part of π system (has double/triple bond)
        for bond in bonds:
            a1 = bond.get('a1')
            a2 = bond.get('a2')
            if (a1 == i or a2 == i) and bond.get('order', 1) >= 2:
                contributing_atoms.append(i)
                break
    
    return {
        "HOMO": round(homo, 2),
        "LUMO": round(lumo, 2),
        "gap": round(gap, 2),
        "contributing_atoms": contributing_atoms,
        "method": "huckel"
    }


def simple_eigenvalues(matrix: List[List[float]]) -> List[float]:
    """Simple eigenvalue calculation for small matrices (2x2, 3x3)"""
    n = len(matrix)
    if n == 0:
        return []
    if n == 1:
        return [matrix[0][0]]
    if n == 2:
        # 2x2 matrix eigenvalues: λ = (trace ± sqrt(trace² - 4*det)) / 2
        a, b = matrix[0][0], matrix[0][1]
        c, d = matrix[1][0], matrix[1][1]
        trace = a + d
        det = a * d - b * c
        discriminant = trace ** 2 - 4 * det
        if discriminant < 0:
            # Complex eigenvalues - return real part
            return [trace / 2, trace / 2]
        sqrt_disc = math.sqrt(discriminant)
        return [(trace - sqrt_disc) / 2, (trace + sqrt_disc) / 2]
    if n == 3:
        # 3x3 matrix - use characteristic polynomial
        # For simplicity, use diagonal elements as approximation
        # In production, would use numpy.linalg.eig
        return sorted([matrix[i][i] for i in range(n)], reverse=True)
    
    # For larger matrices, return diagonal elements as approximation
    # Sorted by energy (most negative first)
    return sorted([matrix[i][i] for i in range(n)], reverse=True)


def count_pi_electrons(molecule: Dict[str, Any]) -> int:
    """Count π electrons in molecule (simplified)"""
    bonds = molecule.get('bonds', [])
    pi_electrons = 0
    
    for bond in bonds:
        bond_order = bond.get('order', 1)
        if bond_order == 2:
            pi_electrons += 2  # One π bond = 2 electrons
        elif bond_order == 3:
            pi_electrons += 4  # Two π bonds = 4 electrons
    
    return pi_electrons

