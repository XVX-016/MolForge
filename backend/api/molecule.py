"""
Molecule operations API endpoints

Phase 6: RDKit Backend Integration
Phase 8: 2D Layout Generation

Provides endpoints for:
- SMILES generation
- MolBlock generation
- Hydrogen normalization
- RDKit validation
- 2D coordinate generation
"""

from typing import List, Dict, Any, Optional
import uuid
import logging
from rdkit import Chem
from rdkit.Chem import rdDistGeom
from backend.chemistry.rdkit_props import compute_properties
from backend.chemistry.errors import ChemistryError
from backend.chemistry.optimize import analyze_structure, suggest_optimizations
from backend.chemistry.diff import calculate_structural_diff
from backend.services.version_service import VersionService
from backend.core.dependencies import get_db
from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from sqlmodel import Session

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/molecule", tags=["molecule"])


class MoleculeRequest(BaseModel):
    """Request format for molecule operations."""
    molecule: Dict[str, Any]  # {atoms: [...], bonds: [...]}
    canonicalize: Optional[bool] = True


class ValidateRequest(BaseModel):
    """Request for molecule validation."""
    molecule: Dict[str, Any]
    smiles: Optional[str] = None  # Alternative: provide SMILES directly


class LayoutRequest(BaseModel):
    """Request for 2D layout generation."""
    molecule: Dict[str, Any]
    smiles: Optional[str] = None  # Alternative: provide SMILES directly
    method: Optional[str] = "coordgen"  # "coordgen" or "rdkit"
    spacing: Optional[float] = 1.5  # Bond length in Angstroms


class ThreeDRequest(BaseModel):
    """Request for 3D coordinate generation."""
    molecule: Dict[str, Any]
    smiles: Optional[str] = None  # Alternative: provide SMILES directly
    optimize: Optional[bool] = True  # Whether to optimize geometry
    method: Optional[str] = "etkdg"  # "etkdg" or "embed"


class PropertyRequest(BaseModel):
    """Request for molecule properties."""
    smiles: str

class CommitRequest(BaseModel):
    """Request to commit a molecule version."""
    smiles: str
    json_graph: Dict[str, Any]
    parent_version_id: Optional[uuid.UUID] = None
    metadata: Optional[Dict[str, Any]] = None


class DashboardRequest(BaseModel):
    """Request for the Studio Dashboard."""
    baseline_version_id: uuid.UUID
    proposal_version_id: Optional[uuid.UUID] = None


class DashboardPayload(BaseModel):
    """Refined Studio Dashboard Response."""
    baseline: Dict[str, Any]
    proposal: Optional[Dict[str, Any]] = None
    diff: Dict[str, Any]
    alerts: List[Dict[str, Any]]
    property_delta: Dict[str, float]
    radar: Dict[str, Dict[str, float]]
    optimization_context: Dict[str, Any]


@router.post("/to-smiles")
async def to_smiles(request: MoleculeRequest):
    """
    Convert molecule to SMILES string.
    
    Uses RDKit for accurate SMILES generation.
    """
    try:
        mol = molecule_dict_to_rdkit(request.molecule)
        if not mol:
            raise HTTPException(status_code=400, detail="Invalid molecule structure")
        
        if request.canonicalize:
            smiles = Chem.MolToSmiles(mol)
        else:
            smiles = Chem.MolToSmiles(mol, canonical=False)
        
        return {
            "smiles": smiles,
            "canonical": request.canonicalize,
        }
    except Exception as e:
        logger.error(f"Error generating SMILES: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/properties")
async def get_properties(request: PropertyRequest):
    """
    Get molecular properties using RDKit.
    """
    try:
        return compute_properties(request.smiles)
    except ChemistryError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Unexpected error in /properties: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail="Internal chemistry engine error")


@router.post("/commit")
async def commit_version(request: CommitRequest, db: Session = Depends(get_db)):
    """
    Commit a molecular draft to immutable persistence.
    """
    try:
        version = VersionService.create_version(
            db=db,
            smiles=request.smiles,
            json_graph=request.json_graph,
            parent_version_id=request.parent_version_id,
            additional_metadata=request.metadata
        )
        return {
            "version_id": version.id,
            "molecule_id": version.molecule_id,
            "version_index": version.version_index,
            "canonical_smiles": version.canonical_smiles,
            "properties": version.properties
        }
    except ChemistryError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Commit Error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail="Failed to commit molecular version")


@router.post("/analyze")
async def analyze_molecule(request: PropertyRequest):
    """
    Perform medicinal chemistry analysis and suggest optimizations.
    """
    try:
        mol = Chem.MolFromSmiles(request.smiles)
        if not mol:
            raise HTTPException(status_code=400, detail="Invalid SMILES")
            
        issues, radar = analyze_structure(mol)
        suggestions = suggest_optimizations(mol)
        
        return {
            "issues": issues,
            "suggestions": suggestions,
            "radar": radar
        }
    except Exception as e:
        logger.error(f"Analysis Error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail="Failed to analyze structure")


@router.get("/diff")
async def get_molecule_diff(base_id: str, prop_id: str, db: Session = Depends(get_db)):
    """
    Compute structural diff between two molecular versions.
    """
    try:
        base_version = db.get(MoleculeVersion, base_id)
        prop_version = db.get(MoleculeVersion, prop_id)
        
        if not base_version or not prop_version:
            raise HTTPException(status_code=404, detail="One or both versions not found")
            
        mol_a = Chem.MolFromSmiles(base_version.canonical_smiles)
        mol_b = Chem.MolFromSmiles(prop_version.canonical_smiles)
        
        if not mol_a or not mol_b:
            raise HTTPException(status_code=400, detail="Could not parse SMILES from versions")
            
        diff = calculate_structural_diff(mol_a, mol_b)
        return diff
    except Exception as e:
        logger.error(f"Diff Endpoint Error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail="Failed to compute structural diff")


@router.post("/dashboard", response_model=DashboardPayload)
async def post_molecule_dashboard(request: DashboardRequest, db: Session = Depends(get_db)):
    """
    Authoritative Studio Dashboard endpoint.
    Orchestrates RDKit properties, structural diffs, and optimization rules.
    """
    try:
        from backend.chemistry.models import MoleculeVersion
        
        base = db.get(MoleculeVersion, request.baseline_version_id)
        if not base:
            raise HTTPException(status_code=404, detail="Baseline version not found")
            
        opt = None
        if request.proposal_version_id:
            opt = db.get(MoleculeVersion, request.proposal_version_id)
            if not opt:
                raise HTTPException(status_code=404, detail="Proposal version not found")
        
        mol_base = Chem.MolFromSmiles(base.canonical_smiles)
        mol_opt = Chem.MolFromSmiles(opt.canonical_smiles) if opt else None
        
        # 1. Properties & Analysis
        # Base properties are usually stored, but we can recompute if needed for consistency
        base_props = base.properties
        
        # Analyze the proposal (or baseline if no proposal)
        analysis_target = mol_opt if mol_opt else mol_base
        issues, radar_raw = analyze_structure(analysis_target)
        suggestions = suggest_optimizations(analysis_target)
        
        # 2. Diffing (if proposal exists)
        diff_payload = {"atoms": {"added": [], "removed": [], "modified": []}, "bonds": {"added": [], "removed": []}}
        prop_delta = {}
        radar_payload = {"baseline": radar_raw, "proposal": {}}
        
        if mol_opt and opt:
            # Structural Diff
            raw_diff = calculate_structural_diff(mol_base, mol_opt)
            # Transform to spec format
            for item in raw_diff.get("proposal", {}).get("atoms", []):
                if item["status"] == "added":
                    diff_payload["atoms"]["added"].append(item["index"])
            for item in raw_diff.get("baseline", {}).get("atoms", []):
                if item["status"] == "deleted":
                    diff_payload["atoms"]["removed"].append(item["index"])
            
            for item in raw_diff.get("proposal", {}).get("bonds", []):
                if item["status"] == "added":
                    # Bond indices are trickier; in RDKit we use atom pairs often
                    # For now, store indices of atoms for the bond
                    bond = mol_opt.GetBondWithIdx(item["index"])
                    diff_payload["bonds"]["added"].append([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
            for item in raw_diff.get("baseline", {}).get("bonds", []):
                if item["status"] == "deleted":
                    bond = mol_base.GetBondWithIdx(item["index"])
                    diff_payload["bonds"]["removed"].append([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])

            # Property Delta
            opt_props = opt.properties
            relevant_keys = ["logp", "tpsa", "qed", "molecular_weight"]
            for key in relevant_keys:
                base_val = base_props.get(key, 0)
                opt_val = opt_props.get(key, 0)
                prop_delta[key] = round(opt_val - base_val, 3)
            
            # Radar Comparison
            _, opt_radar = analyze_structure(mol_opt)
            _, base_radar = analyze_structure(mol_base)
            radar_payload = {
                "baseline": base_radar,
                "proposal": opt_radar
            }

        return {
            "baseline": {
                "version_id": str(base.id),
                "smiles": base.canonical_smiles,
                "properties": base_props
            },
            "proposal": {
                "version_id": str(opt.id),
                "smiles": opt.canonical_smiles,
                "properties": opt.properties
            } if opt else None,
            "diff": diff_payload,
            "alerts": issues,
            "property_delta": prop_delta,
            "radar": radar_payload,
            "optimization_context": {
                "available_rules": suggestions
            }
        }
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Authoritative Dashboard Error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail="Failed to load studio dashboard data")


@router.post("/to-molblock")
async def to_molblock(request: MoleculeRequest):
    """
    Convert molecule to MolBlock format (V2000).
    """
    try:
        mol = molecule_dict_to_rdkit(request.molecule)
        if not mol:
            raise HTTPException(status_code=400, detail="Invalid molecule structure")
        
        molblock = Chem.MolToMolBlock(mol)
        
        return {
            "molblock": molblock,
        }
    except Exception as e:
        logger.error(f"Error generating MolBlock: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/normalize-hydrogens")
async def normalize_hydrogens(request: MoleculeRequest):
    """
    Normalize hydrogens in molecule (add implicit hydrogens).
    """
    try:
        mol = molecule_dict_to_rdkit(request.molecule)
        if not mol:
            raise HTTPException(status_code=400, detail="Invalid molecule structure")
        
        # Add hydrogens
        mol = Chem.AddHs(mol)
        
        # Convert back to molecule dict
        molecule_dict = rdkit_to_molecule_dict(mol)
        
        return {
            "molecule": molecule_dict,
        }
    except Exception as e:
        logger.error(f"Error normalizing hydrogens: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/generate-2d-layout")
async def generate_2d_layout(request: LayoutRequest):
    """
    Generate 2D coordinates for molecule.
    
    Uses RDKit's coordgen or standard 2D coordinate generation.
    Returns molecule with updated atom positions.
    """
    try:
        # Get molecule from dict or SMILES
        mol = None
        if request.smiles:
            mol = Chem.MolFromSmiles(request.smiles)
            if not mol:
                raise HTTPException(status_code=400, detail="Invalid SMILES string")
        else:
            mol = molecule_dict_to_rdkit(request.molecule)
            if not mol:
                raise HTTPException(status_code=400, detail="Invalid molecule structure")
        
        # Generate 2D coordinates
        if request.method == "coordgen":
            # Use CoordGen (better for complex molecules)
            try:
                rdDepictor.Compute2DCoords(mol)
            except:
                # Fallback to standard method
                AllChem.Compute2DCoords(mol)
        else:
            # Use standard RDKit 2D coordinate generation
            AllChem.Compute2DCoords(mol)
        
        # Optimize layout
        try:
            # Try to improve layout with additional optimization
            rdDepictor.Compute2DCoords(mol, clearConfs=True)
        except:
            pass
        
        # Convert back to molecule dict with 2D coordinates
        molecule_dict = rdkit_to_molecule_dict(mol)
        
        return {
            "molecule": molecule_dict,
            "method": request.method,
        }
    except Exception as e:
        logger.error(f"Error generating 2D layout: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/generate-3d")
async def generate_3d_coordinates(request: ThreeDRequest):
    """
    Generate 3D coordinates for molecule.
    
    Uses RDKit ETKDG (Experimental-Torsion-angle preference with Distance Geometry)
    or standard embedding, with optional geometry optimization.
    Returns molecule with updated 3D atom positions.
    """
    try:
        # Get molecule from dict or SMILES
        mol = None
        if request.smiles:
            mol = Chem.MolFromSmiles(request.smiles)
            if not mol:
                raise HTTPException(status_code=400, detail="Invalid SMILES string")
        else:
            mol = molecule_dict_to_rdkit(request.molecule)
            if not mol:
                raise HTTPException(status_code=400, detail="Invalid molecule structure")
        
        # Add hydrogens for 3D
        mol = Chem.AddHs(mol)
        
        # Generate 3D coordinates
        if request.method == "etkdg":
            # Use ETKDG (better for complex molecules)
            try:
                params = rdDistGeom.ETKDGv3()
                params.randomSeed = 42
                embed_result = rdDistGeom.EmbedMolecule(mol, params)
                if embed_result == -1:
                    # Fallback to standard embedding
                    embed_result = AllChem.EmbedMolecule(mol, randomSeed=42)
            except:
                # Fallback to standard embedding
                embed_result = AllChem.EmbedMolecule(mol, randomSeed=42)
        else:
            # Use standard embedding
            embed_result = AllChem.EmbedMolecule(mol, randomSeed=42)
        
        if embed_result == -1:
            # Try with random coords
            embed_result = AllChem.EmbedMolecule(mol, useRandomCoords=True)
            if embed_result == -1:
                raise HTTPException(
                    status_code=500,
                    detail="Failed to generate 3D coordinates for molecule"
                )
        
        # Optimize geometry if requested
        if request.optimize:
            try:
                # Try MMFF first (more accurate)
                AllChem.MMFFOptimizeMolecule(mol)
            except:
                try:
                    # Fallback to UFF (universal)
                    AllChem.UFFOptimizeMolecule(mol)
                except:
                    # If both fail, continue with unoptimized coordinates
                    pass
        
        # Convert back to molecule dict with 3D coordinates
        molecule_dict = rdkit_to_molecule_dict(mol)
        
        return {
            "molecule": molecule_dict,
            "method": request.method,
            "optimized": request.optimize,
        }
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error generating 3D coordinates: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/from-smiles")
async def from_smiles(request: Dict[str, Any]):
    """
    Load molecule from SMILES string.
    
    Returns molecule dict with 2D coordinates.
    """
    try:
        smiles = request.get("smiles")
        if not smiles:
            raise HTTPException(status_code=400, detail="SMILES string required")
        
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise HTTPException(status_code=400, detail="Invalid SMILES string")
        
        # Generate 2D coordinates
        try:
            rdDepictor.Compute2DCoords(mol)
        except:
            AllChem.Compute2DCoords(mol)
        
        # Convert to molecule dict
        molecule_dict = rdkit_to_molecule_dict(mol)
        
        return {
            "molecule": molecule_dict,
        }
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error loading from SMILES: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/from-molblock")
async def from_molblock(request: Dict[str, Any]):
    """
    Load molecule from MolBlock string.
    
    Returns molecule dict with coordinates from MolBlock.
    """
    try:
        molblock = request.get("molblock")
        if not molblock:
            raise HTTPException(status_code=400, detail="MolBlock string required")
        
        mol = Chem.MolFromMolBlock(molblock)
        if not mol:
            raise HTTPException(status_code=400, detail="Invalid MolBlock string")
        
        # Convert to molecule dict (preserves coordinates from MolBlock)
        molecule_dict = rdkit_to_molecule_dict(mol)
        
        return {
            "molecule": molecule_dict,
        }
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error loading from MolBlock: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/validate")
async def validate_molecule(request: ValidateRequest):
    """
    Validate molecule with RDKit.
    
    Returns:
    - valid: bool
    - sanitized_smiles: str (if valid)
    - molblock: str (if valid)
    - errors: List[str]
    """
    try:
        errors: List[str] = []
        
        # Try to create molecule from dict or SMILES
        mol = None
        if request.smiles:
            mol = Chem.MolFromSmiles(request.smiles)
            if not mol:
                errors.append(f"Invalid SMILES: {request.smiles}")
        else:
            mol = molecule_dict_to_rdkit(request.molecule)
            if not mol:
                errors.append("Invalid molecule structure")
        
        if not mol:
            return {
                "valid": False,
                "errors": errors,
            }
        
        # Try to sanitize and canonicalize
        from backend.chemistry.canonical import canonicalize
        try:
            identity = canonicalize(request.smiles if request.smiles else Chem.MolToSmiles(mol))
            
            return {
                "valid": True,
                "sanitized_smiles": identity["canonical_smiles"],
                "inchikey": identity["inchikey"],
                "inchi": identity["inchi"],
                "molblock": Chem.MolToMolBlock(identity["mol"]),
                "errors": [],
            }
        except Exception as sanitize_error:
            errors.append(f"Validation failed: {str(sanitize_error)}")
            return {
                "valid": False,
                "errors": errors,
            }
            
    except Exception as e:
        logger.error(f"Error validating molecule: {e}", exc_info=True)
        return {
            "valid": False,
            "errors": [str(e)],
        }


def molecule_dict_to_rdkit(molecule_dict: Dict[str, Any]) -> Optional[Chem.Mol]:
    """
    Convert molecule dict to RDKit Mol object.
    
    Expected format:
    {
        "atoms": [{"id": "...", "element": "C", "position": [x, y, z], ...}, ...],
        "bonds": [{"id": "...", "atom1": "...", "atom2": "...", "order": 1}, ...]
    }
    """
    try:
        atoms = molecule_dict.get("atoms", [])
        bonds = molecule_dict.get("bonds", [])
        
        if not atoms:
            return None
        
        # Create RDKit molecule
        mol = Chem.RWMol()
        
        # Map atom IDs to RDKit atom indices
        atom_id_to_idx: Dict[str, int] = {}
        
        # Add atoms
        for atom_data in atoms:
            element = atom_data.get("element", "C")
            # Get atomic number from element symbol
            try:
                atomic_num = Chem.GetPeriodicTable().GetAtomicNumber(element)
            except:
                # Fallback: use common elements
                element_map = {
                    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
                    "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
                    "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Br": 35, "I": 53,
                }
                atomic_num = element_map.get(element, 6)  # Default to Carbon
            
            rdkit_atom = Chem.Atom(atomic_num)
            
            # Set charge if provided
            if "charge" in atom_data or "formalCharge" in atom_data:
                charge = atom_data.get("formalCharge") or atom_data.get("charge") or 0
                rdkit_atom.SetFormalCharge(int(charge))
            
            idx = mol.AddAtom(rdkit_atom)
            atom_id_to_idx[atom_data["id"]] = idx
        
        # Add bonds
        for bond_data in bonds:
            atom1_id = bond_data.get("atom1")
            atom2_id = bond_data.get("atom2")
            order = bond_data.get("order", 1)
            
            if atom1_id not in atom_id_to_idx or atom2_id not in atom_id_to_idx:
                continue
            
            idx1 = atom_id_to_idx[atom1_id]
            idx2 = atom_id_to_idx[atom2_id]
            
            # Convert bond order
            if order == 1.5:
                # Aromatic bond
                mol.AddBond(idx1, idx2, Chem.BondType.AROMATIC)
            elif order == 1:
                mol.AddBond(idx1, idx2, Chem.BondType.SINGLE)
            elif order == 2:
                mol.AddBond(idx1, idx2, Chem.BondType.DOUBLE)
            elif order == 3:
                mol.AddBond(idx1, idx2, Chem.BondType.TRIPLE)
            else:
                mol.AddBond(idx1, idx2, Chem.BondType.SINGLE)
        
        # Convert to regular Mol
        mol = mol.GetMol()
        
        # Try to sanitize
        try:
            Chem.SanitizeMol(mol)
        except:
            # If sanitization fails, return unsanitized mol
            # (caller can handle errors)
            pass
        
        return mol
        
    except Exception as e:
        logger.error(f"Error converting molecule dict to RDKit: {e}", exc_info=True)
        return None


def rdkit_to_molecule_dict(mol: Chem.Mol) -> Dict[str, Any]:
    """
    Convert RDKit Mol to molecule dict format.
    Extracts 2D or 3D coordinates from conformer.
    """
    atoms = []
    bonds = []
    
    # Get conformer (prefer 2D, fallback to 3D or default)
    conf = None
    if mol.GetNumConformers() > 0:
        # Prefer conformer 0 (usually 2D if generated)
        conf = mol.GetConformer(0)
    
    # Add atoms
    for i, atom in enumerate(mol.GetAtoms()):
        element = atom.GetSymbol()
        
        # Get position from conformer or default to [0, 0, 0]
        if conf:
            pos = conf.GetAtomPosition(i)
            position = [pos.x, pos.y, pos.z if pos.z else 0.0]
        else:
            position = [0.0, 0.0, 0.0]
        
        atoms.append({
            "id": f"atom_{i}",
            "element": element,
            "position": position,
            "charge": atom.GetFormalCharge(),
            "formalCharge": atom.GetFormalCharge(),
        })
    
    # Add bonds
    for i, bond in enumerate(mol.GetBonds()):
        order = bond.GetBondTypeAsDouble()
        bonds.append({
            "id": f"bond_{i}",
            "atom1": f"atom_{bond.GetBeginAtomIdx()}",
            "atom2": f"atom_{bond.GetEndAtomIdx()}",
            "order": int(order) if order == int(order) else 1.5,
        })
    
    return {
        "atoms": atoms,
        "bonds": bonds,
    }
