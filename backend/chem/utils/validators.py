"""
SMILES validation utilities
"""

try:
    from backend.chemistry.rdkit_props import validate_smiles as _validate_smiles
except ImportError:
    # Fallback implementation
    def _validate_smiles(smiles: str) -> bool:
        """Fallback SMILES validation using RDKit directly."""
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        except ImportError:
            # If RDKit not available, basic string check
            return bool(smiles and len(smiles.strip()) > 0)

# Re-export for convenience
def validate_smiles(smiles: str) -> bool:
    """
    Validate SMILES string.
    
    Args:
        smiles: SMILES string to validate
    
    Returns:
        True if valid, False otherwise
    """
    return _validate_smiles(smiles)

