from rdkit import Chem
from rdkit.Chem import rdFMCS
from typing import Dict, List, Any, Optional
import logging

logger = logging.getLogger(__name__)

def calculate_structural_diff(mol_a: Chem.Mol, mol_b: Chem.Mol) -> Dict[str, Any]:
    """
    Computes a structural diff between two molecules using Maximum Common Substructure (MCS).
    Returns mapping of atoms and bonds to categories: unchanged, added, deleted, modified.
    
    mol_a: Baseline molecule
    mol_b: Proposal/Versioned molecule
    """
    try:
        # 1. Compute MCS
        # We compare elements and bond orders for a high-fidelity diff
        res = rdFMCS.FindMCS([mol_a, mol_b], 
                             atomCompare=rdFMCS.AtomCompare.CompareElements, 
                             bondCompare=rdFMCS.BondCompare.CompareOrder,
                             ringMatchesRingOnly=True)
        
        mcs_mol = Chem.MolFromSmarts(res.smartsString)
        
        # 2. Get mappings from molecules to MCS
        match_a = mol_a.GetSubstructMatch(mcs_mol)
        match_b = mol_b.GetSubstructMatch(mcs_mol)
        
        # 3. Categorize Atoms
        atom_diff_a = [] # Baseline view (highlights deletions)
        atom_diff_b = [] # Proposal view (highlights additions)
        
        # Mapping helpers
        mcs_to_a = {mcs_idx: a_idx for mcs_idx, a_idx in enumerate(match_a)}
        mcs_to_b = {mcs_idx: b_idx for mcs_idx, b_idx in enumerate(match_b)}
        a_to_mcs = {a_idx: mcs_idx for mcs_idx, a_idx in enumerate(match_a)}
        b_to_mcs = {b_idx: mcs_idx for mcs_idx, b_idx in enumerate(match_b)}

        # Baseline Perspective
        for atom in mol_a.GetAtoms():
            idx = atom.GetIdx()
            if idx in a_to_mcs:
                atom_diff_a.append({"index": idx, "status": "unchanged"})
            else:
                atom_diff_a.append({"index": idx, "status": "deleted"})
                
        # Proposal Perspective
        for atom in mol_b.GetAtoms():
            idx = atom.GetIdx()
            if idx in b_to_mcs:
                atom_diff_b.append({"index": idx, "status": "unchanged"})
            else:
                atom_diff_b.append({"index": idx, "status": "added"})

        # 4. Categorize Bonds
        bond_diff_a = []
        bond_diff_b = []

        # Helper to check if bond is in MCS
        def is_bond_in_mcs(bond: Chem.Bond, atom_mapping: Dict[int, int]) -> bool:
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in atom_mapping and a2 in atom_mapping:
                # Potential overlap, check if MCS molecule has bond between these MCS-indices
                m1, m2 = atom_mapping[a1], atom_mapping[a2]
                return mcs_mol.GetBondBetweenAtoms(m1, m2) is not None
            return False

        for bond in mol_a.GetBonds():
            if is_bond_in_mcs(bond, a_to_mcs):
                bond_diff_a.append({"index": bond.GetIdx(), "status": "unchanged"})
            else:
                bond_diff_a.append({"index": bond.GetIdx(), "status": "deleted"})

        for bond in mol_b.GetBonds():
            if is_bond_in_mcs(bond, b_to_mcs):
                bond_diff_b.append({"index": bond.GetIdx(), "status": "unchanged"})
            else:
                bond_diff_b.append({"index": bond.GetIdx(), "status": "added"})

        return {
            "baseline": {
                "atoms": atom_diff_a,
                "bonds": bond_diff_a
            },
            "proposal": {
                "atoms": atom_diff_b,
                "bonds": bond_diff_b
            },
            "mcs_smarts": res.smartsString
        }

    except Exception as e:
        logger.error(f"Structural Diff Error: {e}", exc_info=True)
        return {"error": str(e)}
