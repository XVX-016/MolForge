from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, QED
from backend.chemistry.errors import PropertyComputationError
from backend.chemistry.canonical import canonicalize
from backend.chemistry.cache import descriptor_cache
import logging

logger = logging.getLogger(__name__)

def compute_properties(smiles: str) -> dict:
    """
    Computes molecular properties using RDKit.
    This is the DETERMINISTIC source of truth for MolForge.
    Enforces canonicalization and utilizes caching.
    """
    # 1. Canonicalize at the boundary
    # This gives us the Mol object and the InChIKey for caching
    identity = canonicalize(smiles)
    inchikey = identity["inchikey"]
    
    # 2. Check Cache
    cached = descriptor_cache.get(inchikey)
    if cached:
        logger.debug(f"Cache hit for {inchikey}")
        return cached

    # 3. Compute (if not cached)
    try:
        mol = identity["mol"]
        
        props = {
            # Identity Headers
            "canonical_smiles": identity["canonical_smiles"],
            "inchikey": identity["inchikey"],
            "inchi": identity["inchi"],
            
            # Physics/Chemistry Descriptors
            "molecular_weight": round(Descriptors.MolWt(mol), 3),
            "logp": round(Descriptors.MolLogP(mol), 3),
            "tpsa": round(Descriptors.TPSA(mol), 3),
            "hbd": Lipinski.NumHDonors(mol),
            "hba": Lipinski.NumHAcceptors(mol),
            "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "heavy_atom_count": mol.GetNumHeavyAtoms(),
            "qed": round(QED.qed(mol), 3),
            "complexity": round(Descriptors.BertzCT(mol), 1), # Using Bertz complexity
            "lipinski_violations": int(
                (Descriptors.MolWt(mol) > 500) +
                (Descriptors.MolLogP(mol) > 5) +
                (Lipinski.NumHDonors(mol) > 5) +
                (Lipinski.NumHAcceptors(mol) > 10)
            )
        }
        
        # 4. Burn to Cache
        descriptor_cache.set(inchikey, props)
        
        return props
        
    except Exception as e:
        logger.error(f"RDKit Computation Error for {smiles}: {e}")
        raise PropertyComputationError(f"Could not compute properties: {str(e)}")

def validate_smiles(smiles: str) -> bool:
    """Checks if a SMILES string is valid using RDKit."""
    try:
        # We use canonicalize here to ensure it passes our internal rigorous check
        canonicalize(smiles)
        return True
    except:
        return False
        
def get_canonical_form(smiles: str) -> str:
    """Return only the canonical SMILES string."""
    return canonicalize(smiles)["canonical_smiles"]
