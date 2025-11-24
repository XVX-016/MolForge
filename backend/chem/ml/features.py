"""
Molecular Feature Extraction
Computes fingerprints and descriptors for ML models

Feature Categories:
1. Graph-based features (node/edge attributes)
2. Chemical descriptors (MW, LogP, TPSA, etc.)
3. Engine-derived features (IR, NMR, Energy, Quantum)
4. Optional embeddings (pre-trained molecular representations)
"""
from typing import Dict, List, Any, Optional
from ..kab.utils import get_functional_groups, calculate_surface_area


def extract_features(
    molecule: Dict[str, Any],
    include_engine_features: bool = False,
    engine_outputs: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    Extract comprehensive molecular features for ML prediction
    
    Feature Categories:
    1. Graph-based: atom types, bond types, rings, functional groups
    2. Chemical descriptors: MW, LogP, TPSA, H-bond counts
    3. Engine-derived: IR peaks, NMR shifts, energy, quantum (optional)
    
    Args:
        molecule: Dict with 'atoms' and 'bonds' keys
        include_engine_features: Whether to include engine-derived features
        engine_outputs: Optional dict with IR, NMR, Energy, Quantum outputs
        
    Returns:
        {
            "fingerprint": [0, 1, 0, 1, ...],  # Binary fingerprint
            "descriptors": {
                "molecular_weight": 180.16,
                "logp": 2.3,
                "hbd": 2,  # H-bond donors
                "hba": 4,  # H-bond acceptors
                "rotatable_bonds": 3,
                "aromatic_rings": 1,
                "tpsa": 45.0  # Topological Polar Surface Area
            },
            "graph_features": {
                "node_features": [...],  # Per-atom features
                "edge_features": [...]   # Per-bond features
            },
            "engine_features": {...}  # Optional engine-derived features
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
    
    # Graph-based features
    graph_features = extract_graph_features(molecule)
    
    # Calculate TPSA (Topological Polar Surface Area)
    tpsa = calculate_tpsa(molecule)
    
    # Base feature dict
    features = {
        "fingerprint": fingerprint,
        "descriptors": {
            "molecular_weight": round(molecular_weight, 2),
            "logp": round(logp, 2),
            "hbd": hbd,
            "hba": hba,
            "rotatable_bonds": rotatable_bonds,
            "aromatic_rings": aromatic_rings,
            "tpsa": round(tpsa, 2),
            "num_atoms": len(atoms),
            "num_bonds": len(bonds)
        },
        "graph_features": graph_features
    }
    
    # Add engine-derived features if requested
    if include_engine_features and engine_outputs:
        engine_features = extract_engine_features(molecule, engine_outputs)
        features["engine_features"] = engine_features
    
    return features


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


def extract_graph_features(molecule: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extract graph-based features (node and edge attributes)
    
    Returns:
        {
            "node_features": [
                {"atom_idx": 0, "element": "C", "hybridization": "sp3", "charge": 0, "aromatic": False},
                ...
            ],
            "edge_features": [
                {"bond_idx": 0, "a1": 0, "a2": 1, "order": 1, "conjugated": False, "in_ring": False},
                ...
            ]
        }
    """
    atoms = molecule.get('atoms', [])
    bonds = molecule.get('bonds', [])
    
    # Build adjacency list for ring detection
    adj = {i: [] for i in range(len(atoms))}
    for bond in bonds:
        a1 = bond.get('a1')
        a2 = bond.get('a2')
        if a1 is not None and a2 is not None:
            adj[a1].append(a2)
            adj[a2].append(a1)
    
    # Node features
    node_features = []
    for i, atom in enumerate(atoms):
        element = atom.get('element', 'C')
        neighbors = adj.get(i, [])
        
        # Estimate hybridization (simplified)
        hybridization = estimate_hybridization(element, len(neighbors), bonds, i)
        
        # Check aromaticity (simplified - would need proper ring detection)
        aromatic = is_aromatic_atom(molecule, i)
        
        node_features.append({
            "atom_idx": i,
            "element": element,
            "hybridization": hybridization,
            "charge": atom.get('charge', 0),
            "aromatic": aromatic,
            "degree": len(neighbors),
            "electronegativity": get_electronegativity(element)
        })
    
    # Edge features
    edge_features = []
    for idx, bond in enumerate(bonds):
        a1 = bond.get('a1')
        a2 = bond.get('a2')
        order = bond.get('order', 1)
        
        # Check if conjugated (simplified)
        conjugated = is_conjugated_bond(molecule, a1, a2, order)
        
        # Check if in ring (simplified)
        in_ring = is_bond_in_ring(molecule, a1, a2)
        
        edge_features.append({
            "bond_idx": idx,
            "a1": a1,
            "a2": a2,
            "order": order,
            "conjugated": conjugated,
            "in_ring": in_ring
        })
    
    return {
        "node_features": node_features,
        "edge_features": edge_features
    }


def estimate_hybridization(element: str, degree: int, bonds: List[Dict], atom_idx: int) -> str:
    """Estimate atom hybridization"""
    # Count double/triple bonds
    double_bonds = sum(1 for b in bonds 
                      if (b.get('a1') == atom_idx or b.get('a2') == atom_idx) 
                      and b.get('order', 1) >= 2)
    
    if double_bonds > 0:
        return 'sp2' if degree <= 3 else 'sp'
    return 'sp3'


def is_aromatic_atom(molecule: Dict[str, Any], atom_idx: int) -> bool:
    """Simplified aromaticity detection"""
    functional_groups = get_functional_groups(molecule)
    for group in functional_groups:
        if atom_idx in group.get('atoms', []) and group.get('type') == 'aromatic':
            return True
    return False


def get_electronegativity(element: str) -> float:
    """Get Pauling electronegativity"""
    en = {
        'H': 2.1, 'C': 2.5, 'N': 3.0, 'O': 3.5,
        'F': 4.0, 'P': 2.1, 'S': 2.5, 'Cl': 3.0,
        'Br': 2.8, 'I': 2.5
    }
    return en.get(element, 2.5)


def is_conjugated_bond(molecule: Dict[str, Any], a1: int, a2: int, order: int) -> bool:
    """Check if bond is conjugated (simplified)"""
    return order >= 2


def is_bond_in_ring(molecule: Dict[str, Any], a1: int, a2: int) -> bool:
    """Simplified ring detection (would need proper cycle detection)"""
    # Placeholder - would implement proper cycle detection
    return False


def calculate_tpsa(molecule: Dict[str, Any]) -> float:
    """
    Calculate Topological Polar Surface Area (TPSA)
    Simplified version based on functional groups
    """
    functional_groups = get_functional_groups(molecule)
    
    # TPSA contributions (approximate, in Å²)
    tpsa_contributions = {
        'hydroxyl': 20.23,
        'amino': 26.02,
        'carbonyl': 17.07,
        'carboxyl': 37.30,
        'ester': 26.30,
        'amide': 43.69
    }
    
    tpsa = 0.0
    for group in functional_groups:
        group_type = group.get('type', '')
        if group_type in tpsa_contributions:
            tpsa += tpsa_contributions[group_type]
    
    return tpsa


def extract_engine_features(molecule: Dict[str, Any], engine_outputs: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extract features from engine outputs (IR, NMR, Energy, Quantum)
    
    Args:
        molecule: Original molecule dict
        engine_outputs: Dict with keys like 'ir', 'nmr', 'energy', 'quantum'
        
    Returns:
        {
            "ir_features": {"peak_count": 5, "functional_group_peaks": {...}},
            "nmr_features": {"avg_shift": 3.5, "shift_range": 2.0},
            "energy_features": {"total_energy": -120.5, "components": {...}},
            "quantum_features": {"homo_lumo_gap": 4.2, "esp_range": 0.5}
        }
    """
    features = {}
    
    # IR features
    if 'ir' in engine_outputs:
        ir_data = engine_outputs['ir']
        peaks = ir_data.get('peaks', [])
        features['ir_features'] = {
            'peak_count': len(peaks),
            'functional_group_peaks': {}
        }
        # Count peaks per functional group
        for peak in peaks:
            group = peak.get('functional_group', 'unknown')
            features['ir_features']['functional_group_peaks'][group] = \
                features['ir_features']['functional_group_peaks'].get(group, 0) + 1
    
    # NMR features
    if 'nmr' in engine_outputs:
        nmr_data = engine_outputs['nmr']
        shifts = nmr_data.get('shifts', [])
        if shifts:
            shift_values = [s.get('shift', 0) for s in shifts]
            features['nmr_features'] = {
                'avg_shift': sum(shift_values) / len(shift_values),
                'shift_range': max(shift_values) - min(shift_values) if shift_values else 0.0
            }
    
    # Energy features
    if 'energy' in engine_outputs:
        energy_data = engine_outputs['energy']
        features['energy_features'] = {
            'total_energy': energy_data.get('total_energy', 0.0),
            'bond_energy': energy_data.get('bond_energy', 0.0),
            'angle_energy': energy_data.get('angle_energy', 0.0),
            'torsion_energy': energy_data.get('torsion_energy', 0.0)
        }
    
    # Quantum features
    if 'quantum' in engine_outputs:
        quantum_data = engine_outputs['quantum']
        homo_lumo = quantum_data.get('homo_lumo', {})
        esp = quantum_data.get('esp', {})
        
        features['quantum_features'] = {
            'homo_lumo_gap': homo_lumo.get('gap', 0.0),
            'homo_energy': homo_lumo.get('HOMO', 0.0),
            'lumo_energy': homo_lumo.get('LUMO', 0.0)
        }
        
        if esp and 'esp_values' in esp:
            esp_values = [v.get('esp', 0) for v in esp['esp_values']]
            if esp_values:
                features['quantum_features']['esp_range'] = max(esp_values) - min(esp_values)
                features['quantum_features']['esp_avg'] = sum(esp_values) / len(esp_values)
    
    return features

