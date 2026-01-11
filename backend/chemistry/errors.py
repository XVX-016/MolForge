
class ChemistryError(Exception):
    """Base class for all chemistry-related errors."""
    pass

class InvalidSmilesError(ChemistryError):
    """Raised when a SMILES string is invalid."""
    pass

class PropertyComputationError(ChemistryError):
    """Raised when RDKit fails to compute a property."""
    pass
