"""
RDKit Index - ECFP fingerprint computation

Uses existing featurizer to compute fingerprints.
"""

from typing import Optional, Set
try:
    from backend.chem.screening.fingerprint_index import compute_ecfp_fingerprint
except ImportError:
    # Fallback implementation
    compute_ecfp_fingerprint = None


def compute_ecfp(smiles: str, radius: int = 2, n_bits: int = 2048) -> Optional[Set[int]]:
    """
    Compute ECFP fingerprint for SMILES.
    
    Args:
        smiles: SMILES string
        radius: ECFP radius (default: 2 for ECFP4)
        n_bits: Number of bits in fingerprint
    
    Returns:
        Set of active bit indices, or None if invalid SMILES
    """
    if compute_ecfp_fingerprint:
        return compute_ecfp_fingerprint(smiles, radius=radius, n_bits=n_bits)
    
    # Fallback: direct RDKit implementation
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
        return set(fp.GetOnBits())
    except ImportError:
        return None
