"""
Knowledge-based Rule Engine
Applies chemical knowledge rules and ML-based predictions
"""
from typing import Dict, List, Any
from .utils import get_functional_groups, get_atom_neighbors


def query_knowledge_rules(molecule: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Apply knowledge-based rules to molecule
    
    Args:
        molecule: Molecule dict
        
    Returns:
        List of rules with:
        - rule: Rule identifier
        - description: Human-readable description
        - affected_atoms: List of atom indices
        - confidence: Confidence level (0-1)
        - category: Rule category (e.g., 'reactivity', 'stability', 'activity')
    """
    atoms = molecule.get('atoms', [])
    bonds = molecule.get('bonds', [])
    rules = []
    
    # Build adjacency list
    adj = {i: [] for i in range(len(atoms))}
    for bond in bonds:
        a1 = bond.get('a1')
        a2 = bond.get('a2')
        if a1 is not None and a2 is not None:
            adj[a1].append(a2)
            adj[a2].append(a1)
    
    # Rule 1: Aromatic rings increase stability
    aromatic_rings = detect_aromatic_rings(molecule)
    for ring in aromatic_rings:
        rules.append({
            'rule': 'aromatic_stability',
            'description': 'Aromatic ring detected - contributes to molecular stability',
            'affected_atoms': ring['atoms'],
            'confidence': 0.8,
            'category': 'stability'
        })
    
    # Rule 2: Hydrogen bonding capability
    h_bond_donors = []
    h_bond_acceptors = []
    
    for i, atom in enumerate(atoms):
        element = atom.get('element', '')
        neighbors = adj.get(i, [])
        
        # H-bond donors: O-H, N-H
        if element == 'O' or element == 'N':
            for neighbor_idx in neighbors:
                neighbor = atoms[neighbor_idx] if neighbor_idx < len(atoms) else None
                if neighbor and neighbor.get('element') == 'H':
                    h_bond_donors.append(i)
                    break
        
        # H-bond acceptors: O, N with lone pairs
        if element == 'O' or element == 'N':
            # Check if has lone pairs (fewer than max bonds)
            max_bonds = 2 if element == 'O' else 3
            if len(neighbors) < max_bonds:
                h_bond_acceptors.append(i)
    
    if h_bond_donors:
        rules.append({
            'rule': 'h_bond_donor',
            'description': f'{len(h_bond_donors)} hydrogen bond donor(s) detected - can form H-bonds',
            'affected_atoms': h_bond_donors,
            'confidence': 0.9,
            'category': 'reactivity'
        })
    
    if h_bond_acceptors:
        rules.append({
            'rule': 'h_bond_acceptor',
            'description': f'{len(h_bond_acceptors)} hydrogen bond acceptor(s) detected - can accept H-bonds',
            'affected_atoms': h_bond_acceptors,
            'confidence': 0.9,
            'category': 'reactivity'
        })
    
    # Rule 3: Lipinski's Rule of Five (simplified)
    lipinski_violations = check_lipinski_rules(molecule)
    if lipinski_violations:
        for violation in lipinski_violations:
            rules.append({
                'rule': 'lipinski_violation',
                'description': violation['description'],
                'affected_atoms': [],
                'confidence': 0.7,
                'category': 'drug_likeness'
            })
    
    # Rule 4: Functional group reactivity
    functional_groups = get_functional_groups(molecule)
    for group in functional_groups:
        group_type = group.get('type', '')
        
        reactivity_rules = {
            'carbonyl': {
                'description': 'Carbonyl group - susceptible to nucleophilic attack',
                'confidence': 0.8,
                'category': 'reactivity'
            },
            'hydroxyl': {
                'description': 'Hydroxyl group - can participate in substitution reactions',
                'confidence': 0.7,
                'category': 'reactivity'
            },
            'amino': {
                'description': 'Amino group - can act as nucleophile or base',
                'confidence': 0.8,
                'category': 'reactivity'
            }
        }
        
        if group_type in reactivity_rules:
            rule_info = reactivity_rules[group_type]
            rules.append({
                'rule': f'{group_type}_reactivity',
                'description': rule_info['description'],
                'affected_atoms': group.get('atoms', []),
                'confidence': rule_info['confidence'],
                'category': rule_info['category']
            })
    
    # Rule 5: Chirality detection
    chiral_centers = detect_chiral_centers(molecule)
    if chiral_centers:
        rules.append({
            'rule': 'chiral_centers',
            'description': f'{len(chiral_centers)} chiral center(s) detected - molecule may have stereoisomers',
            'affected_atoms': chiral_centers,
            'confidence': 0.6,
            'category': 'stereochemistry'
        })
    
    return rules


def detect_aromatic_rings(molecule: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Detect aromatic rings (simplified - looks for 6-membered rings with alternating bonds)"""
    atoms = molecule.get('atoms', [])
    bonds = molecule.get('bonds', [])
    
    if len(atoms) < 6:
        return []
    
    # Build adjacency list
    adj = {i: [] for i in range(len(atoms))}
    for bond in bonds:
        a1 = bond.get('a1')
        a2 = bond.get('a2')
        if a1 is not None and a2 is not None:
            adj[a1].append(a2)
            adj[a2].append(a1)
    
    # Simple heuristic: 6-membered rings with carbon atoms
    # In reality, aromaticity requires Hückel's rule and proper electron delocalization
    rings = []
    visited = set()
    
    def find_6_ring(start: int, current: int, path: List[int]):
        if len(path) > 6:
            return
        
        for neighbor in adj.get(current, []):
            if neighbor == start and len(path) == 6:
                # Check if mostly carbons
                ring_atoms = [atoms[i] for i in path if i < len(atoms)]
                carbon_count = sum(1 for a in ring_atoms if a.get('element') == 'C')
                if carbon_count >= 4:  # At least 4 carbons
                    rings.append({
                        'atoms': path + [start],
                        'size': 6
                    })
                return
            elif neighbor not in path and neighbor not in visited:
                find_6_ring(start, neighbor, path + [neighbor])
    
    for start in range(len(atoms)):
        if start not in visited:
            find_6_ring(start, start, [start])
            visited.add(start)
    
    return rings


def check_lipinski_rules(molecule: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Check Lipinski's Rule of Five (simplified)"""
    atoms = molecule.get('atoms', [])
    violations = []
    
    # Count atoms for MW estimation
    atomic_weights = {
        'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999,
        'F': 18.998, 'P': 30.974, 'S': 32.065, 'Cl': 35.453,
        'Br': 79.904, 'I': 126.904
    }
    
    molecular_weight = sum(
        atomic_weights.get(atom.get('element', 'C'), 12.0)
        for atom in atoms
    )
    
    # Rule 1: MW < 500 Da
    if molecular_weight > 500:
        violations.append({
            'description': f'Molecular weight ({molecular_weight:.1f} Da) exceeds 500 Da (Rule of Five violation)'
        })
    
    # Rule 2: H-bond donors <= 5 (simplified: count O-H and N-H)
    h_donors = 0
    for atom in atoms:
        element = atom.get('element', '')
        if element == 'O' or element == 'N':
            # Simplified: assume can be H-bond donor if has H neighbor
            # In reality, need to check actual bonding
            h_donors += 1
    
    if h_donors > 5:
        violations.append({
            'description': f'Too many potential H-bond donors ({h_donors} > 5) (Rule of Five violation)'
        })
    
    # Rule 3: H-bond acceptors <= 10 (simplified: count O and N)
    h_acceptors = sum(1 for atom in atoms if atom.get('element') in ['O', 'N'])
    if h_acceptors > 10:
        violations.append({
            'description': f'Too many H-bond acceptors ({h_acceptors} > 10) (Rule of Five violation)'
        })
    
    return violations


def detect_chiral_centers(molecule: Dict[str, Any]) -> List[int]:
    """Detect chiral centers (simplified - carbon with 4 different substituents)"""
    atoms = molecule.get('atoms', [])
    bonds = molecule.get('bonds', [])
    
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
    
    chiral_centers = []
    
    for i, atom in enumerate(atoms):
        if atom.get('element') != 'C':
            continue
        
        neighbors = adj.get(i, [])
        if len(neighbors) != 4:  # Must have 4 bonds
            continue
        
        # Check if substituents are different (simplified: check element types)
        neighbor_elements = [atoms[n].get('element', '') for n in neighbors if n < len(atoms)]
        
        # If all 4 are different, likely chiral
        if len(set(neighbor_elements)) == 4:
            chiral_centers.append(i)
    
    return chiral_centers


def apply_ml_rules(molecule: Dict[str, Any]) -> Dict[str, Any]:
    """
    Apply ML-based property predictions (placeholder for future ML integration)
    
    For now, returns heuristic-based predictions
    """
    atoms = molecule.get('atoms', [])
    functional_groups = get_functional_groups(molecule)
    
    # Heuristic predictions (can be replaced with actual ML model)
    predictions = {
        'drug_likeness': 0.5,  # Default
        'bioavailability': 0.5,
        'toxicity_risk': 0.3,
        'metabolic_stability': 0.5
    }
    
    # Adjust based on molecular properties
    if len(atoms) < 20:
        predictions['drug_likeness'] += 0.1
    elif len(atoms) > 50:
        predictions['drug_likeness'] -= 0.2
    
    # Polar groups increase bioavailability
    polar_groups = [g for g in functional_groups if g.get('type') in ['hydroxyl', 'amino', 'carbonyl']]
    if polar_groups:
        predictions['bioavailability'] += min(0.2, len(polar_groups) * 0.05)
    
    # Clamp values
    for key in predictions:
        predictions[key] = max(0.0, min(1.0, predictions[key]))
    
    return predictions

