from rdkit import Chem
from backend.chemistry.errors import ChemistryError

def canonicalize(smiles: str) -> dict:
    """
    Rigorously canonicalize a SMILES string.
    Returns the Mol object and standard identifiers (InChI, InChIKey).
    """
    if not smiles or not smiles.strip():
        raise ChemistryError("Empty SMILES string")
        
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ChemistryError(f"Invalid SMILES structure: {smiles}")

    try:
        # Generate canonical isomeric SMILES
        canonical_smiles = Chem.MolToSmiles(
            mol, 
            canonical=True, 
            isomericSmiles=True
        )
        
        # Generate standard identity strings
        inchi = Chem.MolToInchi(mol)
        inchikey = Chem.MolToInchiKey(mol)
        
        return {
            "mol": mol,  # Preserve the object to avoid redundant reparsing
            "canonical_smiles": canonical_smiles,
            "inchi": inchi,
            "inchikey": inchikey
        }
    except Exception as e:
        raise ChemistryError(f"RDKit Error during canonicalization: {str(e)}")
