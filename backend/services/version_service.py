from typing import Optional, Dict, Any, List
import uuid
from sqlmodel import Session, select, func
from backend.chemistry.models import MoleculeIdentity, MoleculeVersion
from backend.chemistry.rdkit_props import compute_properties
from backend.chemistry.canonical import canonicalize
import logging

logger = logging.getLogger(__name__)

class VersionService:
    """
    Manages the lifecycle of molecules and their versions.
    Enforces immutability and identity-first persistence.
    """
    
    @staticmethod
    def create_version(
        db: Session,
        smiles: str,
        json_graph: Dict[str, Any],
        user_id: Optional[uuid.UUID] = None,
        parent_version_id: Optional[uuid.UUID] = None,
        state: str = "VERSIONED",
        metadata: Optional[Dict[str, Any]] = None
    ) -> MoleculeVersion:
        """
        Creates a new immutable version of a molecule.
        If the molecule (identity) doesn't exist, it creates it first.
        """
        # 1. Canonicalize and get identifiers (The Guard)
        identity_data = canonicalize(smiles)
        inchikey = identity_data["inchikey"]
        canonical_smiles = identity_data["canonical_smiles"]
        
        # 2. Get or Create the Molecule Identity
        statement = select(MoleculeIdentity).where(MoleculeIdentity.inchikey == inchikey)
        molecule = db.exec(statement).first()
        
        if not molecule:
            molecule = MoleculeIdentity(
                inchikey=inchikey,
                user_id=user_id
            )
            db.add(molecule)
            db.flush() # Get ID without committing
            logger.info(f"Created new molecule identity: {inchikey}")
            
        # 3. Determine Version Index & Parent Enforcement
        # Fetch last version to enforce parent-linkage and index increment
        last_version_stmt = select(MoleculeVersion).where(MoleculeVersion.molecule_id == molecule.id).order_by(MoleculeVersion.version_index.desc())
        last_version = db.exec(last_version_stmt).first()
        
        if last_version:
            # Enforce parent_version_id after first version
            if not parent_version_id:
                logger.warning(f"No parent version provided for molecule {molecule.id}, attaching to last: {last_version.id}")
                parent_version_id = last_version.id
            version_index = last_version.version_index + 1
        else:
            version_index = 1
        
        # 4. Snapshot Properties...
        # Note: In practice, we'd pass computed props to avoid redundant work, 
        # but for safety we compute or fetch from cache here.
        props = compute_properties(canonical_smiles)
        
        # 5. Create the Version
        version = MoleculeVersion(
            molecule_id=molecule.id,
            version_index=version_index,
            parent_version_id=parent_version_id,
            canonical_smiles=canonical_smiles,
            json_graph=json_graph,
            properties=props,
            state=state,
            additional_metadata=metadata or {},
            user_id=user_id
        )
        
        db.add(version)
        db.commit()
        db.refresh(version)
        
        logger.info(f"Committed new version {version_index} for molecule {inchikey}")
        return version

    @staticmethod
    def get_version(db: Session, version_id: uuid.UUID) -> Optional[MoleculeVersion]:
        """Fetch a specific version by ID."""
        return db.get(MoleculeVersion, version_id)

    @staticmethod
    def get_history(db: Session, molecule_id: uuid.UUID) -> List[MoleculeVersion]:
        """Fetch the full history of a molecule."""
        statement = select(MoleculeVersion).where(MoleculeVersion.molecule_id == molecule_id).order_by(MoleculeVersion.version_index.asc())
        return list(db.exec(statement))

def diff_properties(old_props: Dict[str, Any], new_props: Dict[str, Any]) -> Dict[str, Any]:
    """
    Computes the numeric delta between two property sets.
    """
    diff = {}
    
    # We only diff numeric properties that are relevant for medicinal chemistry
    numeric_keys = [
        "molecular_weight", "logp", "tpsa", "hbd", "hba", 
        "rotatable_bonds", "heavy_atom_count", "lipinski_violations"
    ]
    
    for key in numeric_keys:
        old_val = old_props.get(key, 0)
        new_val = new_props.get(key, 0)
        
        diff[key] = {
            "from": old_val,
            "to": new_val,
            "delta": round(new_val - old_val, 3)
        }
        
    return diff
