from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from backend.chemistry.errors import InvalidSmilesError, PropertyComputationError
import logging

logger = logging.getLogger(__name__)

def compute_properties(smiles: str) -> dict:
    """
    Computes molecular properties using RDKit.
    This is the DETERMINISTIC source of truth for MolForge.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise InvalidSmilesError(f"RDKit: Invalid SMILES string '{smiles}'")

        # Basic Sanitization to ensure valence and aromaticity are correct
        try:
            Chem.SanitizeMol(mol)
        except:
            pass # Continue if possible, or raise if critical property fails

        return {
            "molecular_weight": round(Descriptors.MolWt(mol), 3),
            "logp": round(Descriptors.MolLogP(mol), 3),
            "tpsa": round(Descriptors.TPSA(mol), 3),
            "hbd": Lipinski.NumHDonors(mol),
            "hba": Lipinski.NumHAcceptors(mol),
            "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
            "formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
            "heavy_atom_count": mol.GetNumHeavyAtoms(),
            "lipinski_violations": int(
                (Descriptors.MolWt(mol) > 500) +
                (Descriptors.MolLogP(mol) > 5) +
                (Lipinski.NumHDonors(mol) > 5) +
                (Lipinski.NumHAcceptors(mol) > 10)
            )
        }
    except InvalidSmilesError:
        raise
    except Exception as e:
        logger.error(f"RDKit Computation Error for {smiles}: {e}")
        raise PropertyComputationError(f"Could not compute properties: {str(e)}")

def validate_smiles(smiles: str) -> bool:
    """Checks if a SMILES string is valid using RDKit."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except:
        return False
