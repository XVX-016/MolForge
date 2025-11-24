"""
Electrostatic Potential (ESP) Calculation Engine
Computes approximate ESP based on atomic charges
"""
from typing import Dict, List, Any
from ..kab.utils import calculate_distance


def calculate_esp(molecule: Dict[str, Any]) -> Dict[str, Any]:
    """
    Calculate electrostatic potential at each atom
    
    Uses simplified point charge model
    
    Args:
        molecule: Dict with 'atoms' and 'bonds' keys
        
    Returns:
        {
            "esp_values": [
                {"atom": 0, "esp": -0.32},
                {"atom": 1, "esp": 0.11}
            ],
            "method": "point_charge"
        }
    """
    atoms = molecule.get('atoms', [])
    bonds = molecule.get('bonds', [])
    
    if not atoms:
        return {
            "esp_values": [],
            "method": "point_charge"
        }
    
    # Estimate atomic charges using simple electronegativity rules
    charges = estimate_atomic_charges(molecule)
    
    # Calculate ESP at each atom position
    esp_values = []
    for i, atom in enumerate(atoms):
        atom_pos = atom.get('position', [0, 0, 0])
        esp = 0.0
        
        # Sum contributions from all other atoms
        for j, other_atom in enumerate(atoms):
            if i == j:
                continue
            
            other_pos = other_atom.get('position', [0, 0, 0])
            distance = calculate_distance(atom_pos, other_pos)
            
            if distance > 0.01:  # Avoid division by zero
                # ESP = q / r (in atomic units, simplified)
                esp += charges[j] / distance
        
        esp_values.append({
            "atom": i,
            "esp": round(esp, 3)
        })
    
    return {
        "esp_values": esp_values,
        "method": "point_charge"
    }


def estimate_atomic_charges(molecule: Dict[str, Any]) -> List[float]:
    """
    Estimate atomic charges using simple electronegativity rules
    
    Returns list of charges for each atom
    """
    atoms = molecule.get('atoms', [])
    bonds = molecule.get('bonds', [])
    
    if not atoms:
        return []
    
    # Electronegativity values (Pauling scale, simplified)
    electronegativity = {
        'H': 2.1, 'C': 2.5, 'N': 3.0, 'O': 3.5,
        'F': 4.0, 'P': 2.1, 'S': 2.5, 'Cl': 3.0,
        'Br': 2.8, 'I': 2.5
    }
    
    charges = [0.0] * len(atoms)
    
    # Simple charge estimation based on bonds
    for bond in bonds:
        a1 = bond.get('a1')
        a2 = bond.get('a2')
        if a1 is None or a2 is None:
            continue
        
        if a1 >= len(atoms) or a2 >= len(atoms):
            continue
        
        element1 = atoms[a1].get('element', 'C')
        element2 = atoms[a2].get('element', 'C')
        
        en1 = electronegativity.get(element1, 2.5)
        en2 = electronegativity.get(element2, 2.5)
        
        # More electronegative atom gets partial negative charge
        diff = en2 - en1
        partial_charge = diff * 0.1  # Scale factor
        
        charges[a1] -= partial_charge
        charges[a2] += partial_charge
    
    return charges

