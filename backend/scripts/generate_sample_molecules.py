#!/usr/bin/env python3
"""
Generate sample molecules for the library
Creates molecules with SMILES, formula, molfile, and thumbnails
"""
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import base64
import io

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.services.thumbnail_service import ThumbnailService


# Sample molecules with common compounds
SAMPLE_MOLECULES = [
    {
        "name": "Aspirin",
        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "formula": "C9H8O4",
        "description": "Common pain reliever and anti-inflammatory drug"
    },
    {
        "name": "Caffeine",
        "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "formula": "C8H10N4O2",
        "description": "Stimulant found in coffee and tea"
    },
    {
        "name": "Glucose",
        "smiles": "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O",
        "formula": "C6H12O6",
        "description": "Simple sugar, primary energy source"
    },
    {
        "name": "Ethanol",
        "smiles": "CCO",
        "formula": "C2H6O",
        "description": "Alcohol found in beverages"
    },
    {
        "name": "Paracetamol",
        "smiles": "CC(=O)NC1=CC=C(C=C1)O",
        "formula": "C8H9NO2",
        "description": "Pain reliever and fever reducer"
    },
    {
        "name": "Ibuprofen",
        "smiles": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
        "formula": "C13H18O2",
        "description": "Nonsteroidal anti-inflammatory drug"
    },
    {
        "name": "Penicillin G",
        "smiles": "CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C",
        "formula": "C16H18N2O4S",
        "description": "Antibiotic medication"
    },
    {
        "name": "Morphine",
        "smiles": "CN1CC[C@]23C4=C5C=CC(=C2O1)C1=C3C(=O)C=C(C1)O[C@H]4[C@@H](O)C=C5",
        "formula": "C17H19NO3",
        "description": "Opioid pain medication"
    },
    {
        "name": "Vitamin C",
        "smiles": "C([C@@H]([C@H]1C(=C(C(=O)O1)O)O)O)O",
        "formula": "C6H8O6",
        "description": "Ascorbic acid, essential vitamin"
    },
    {
        "name": "Nicotine",
        "smiles": "CN1CCC[C@H]1c2cccnc2",
        "formula": "C10H14N2",
        "description": "Stimulant found in tobacco"
    },
    {
        "name": "Serotonin",
        "smiles": "C1=CC2=C(C=C1O)C(=CN2)CCN",
        "formula": "C10H12N2O",
        "description": "Neurotransmitter"
    },
    {
        "name": "Dopamine",
        "smiles": "C1=CC(=C(C=C1CCN)O)O",
        "formula": "C8H11NO2",
        "description": "Neurotransmitter"
    },
    {
        "name": "Adrenaline",
        "smiles": "CNC[C@H](Cc1ccc(O)c(O)c1)O",
        "formula": "C9H13NO3",
        "description": "Hormone and neurotransmitter"
    },
    {
        "name": "Testosterone",
        "smiles": "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C",
        "formula": "C19H28O2",
        "description": "Male sex hormone"
    },
    {
        "name": "Cholesterol",
        "smiles": "CC(C)CCCC(C)C1CCC2C1(C(C3C2CCC4C3(CCC(C4)O)C)C)C",
        "formula": "C27H46O",
        "description": "Sterol found in cell membranes"
    },
    {
        "name": "LSD",
        "smiles": "CCN(CC)C(=O)[C@H]1CN([C@@H]2Cc3c[nH]c4c3c(ccc4)C2=C1)C",
        "formula": "C20H25N3O",
        "description": "Psychoactive compound"
    },
    {
        "name": "Cocaine",
        "smiles": "COC(=O)[C@H]1C2CCC(CC2OC(=O)C1=O)C3=CC=CC=C3",
        "formula": "C17H21NO4",
        "description": "Stimulant alkaloid"
    },
    {
        "name": "Quinine",
        "smiles": "COC1=CC2=C(C=CN=C2C=C1)[C@H]3C[C@@H]4CCN3CC4C=C",
        "formula": "C20H24N2O2",
        "description": "Antimalarial medication"
    },
    {
        "name": "Codeine",
        "smiles": "CN1CC[C@]23C4=C5C=CC(=C2O1)C1=C3C(=O)C=C(C1)O[C@H]4[C@@H](O)C=C5",
        "formula": "C18H21NO3",
        "description": "Opioid pain medication"
    },
    {
        "name": "Vanillin",
        "smiles": "COc1cc(C=O)ccc1O",
        "formula": "C8H8O3",
        "description": "Vanilla flavor compound"
    }
]


def generate_molecule_data(mol_data: Dict) -> Dict:
    """
    Generate complete molecule data including molfile and thumbnail
    
    Args:
        mol_data: Dictionary with name, smiles, formula, description
        
    Returns:
        Dictionary with all molecule properties
    """
    smiles = mol_data["smiles"]
    name = mol_data["name"]
    
    print(f"Processing {name}...")
    
    try:
        # Generate 3D molfile
        molfile_3d = ThumbnailService.generate_3d_molfile(smiles=smiles)
        
        # Generate thumbnail
        thumbnail_b64 = ThumbnailService.generate_2d_thumbnail_base64(
            smiles=smiles,
            size=(600, 600)
        )
        
        # Generate properties (mock data for now)
        properties = {
            "stability": 0.75 + (hash(smiles) % 25) / 100,  # Random between 0.75-1.0
            "toxicity": 0.1 + (hash(smiles) % 30) / 100,   # Random between 0.1-0.4
            "solubility": 0.5 + (hash(smiles) % 40) / 100,  # Random between 0.5-0.9
            "bioavailability": 0.6 + (hash(smiles) % 30) / 100,  # Random between 0.6-0.9
            "novelty": 0.3 + (hash(smiles) % 50) / 100     # Random between 0.3-0.8
        }
        
        return {
            "name": name,
            "smiles": smiles,
            "formula": mol_data["formula"],
            "molfile": molfile_3d,
            "thumbnail_b64": thumbnail_b64,
            "properties": json.dumps(properties),
            "description": mol_data.get("description", ""),
            "stability_score": properties["stability"]
        }
        
    except Exception as e:
        print(f"  ⚠️  Error processing {name}: {e}")
        # Return minimal data if generation fails
        return {
            "name": name,
            "smiles": smiles,
            "formula": mol_data["formula"],
            "molfile": None,
            "thumbnail_b64": None,
            "properties": None,
            "description": mol_data.get("description", ""),
            "stability_score": None
        }


def generate_all_molecules() -> List[Dict]:
    """Generate data for all sample molecules"""
    results = []
    
    print(f"Generating {len(SAMPLE_MOLECULES)} sample molecules...\n")
    
    for mol_data in SAMPLE_MOLECULES:
        result = generate_molecule_data(mol_data)
        results.append(result)
        print(f"  ✓ {result['name']} - {result['formula']}")
    
    return results


def save_to_json(output_path: str = "sample_molecules.json"):
    """Save generated molecules to JSON file"""
    molecules = generate_all_molecules()
    
    output_file = Path(output_path)
    output_file.write_text(json.dumps(molecules, indent=2))
    
    print(f"\n✅ Saved {len(molecules)} molecules to {output_path}")
    return molecules


def print_supabase_insert_sql(molecules: List[Dict], user_id: str = "YOUR_USER_ID"):
    """Print SQL INSERT statements for Supabase"""
    print("\n" + "="*60)
    print("Supabase INSERT SQL (replace YOUR_USER_ID with actual user ID):")
    print("="*60 + "\n")
    
    for mol in molecules:
        # Escape single quotes in strings
        name = mol["name"].replace("'", "''")
        smiles = mol["smiles"].replace("'", "''") if mol["smiles"] else "NULL"
        formula = mol["formula"].replace("'", "''") if mol["formula"] else "NULL"
        molfile = mol["molfile"].replace("'", "''") if mol["molfile"] else "NULL"
        thumbnail = mol["thumbnail_b64"].replace("'", "''") if mol["thumbnail_b64"] else "NULL"
        properties = mol["properties"].replace("'", "''") if mol["properties"] else "NULL"
        description = mol.get("description", "").replace("'", "''")
        stability = mol.get("stability_score") or "NULL"
        
        sql = f"""INSERT INTO molecules (name, smiles, formula, molfile, thumbnail_b64, properties, user_id, stability_score, created_at)
VALUES ('{name}', '{smiles}', '{formula}', {'NULL' if not molfile or len(molfile) > 10000 else f"'{molfile}'"}, 
        {'NULL' if not thumbnail or len(thumbnail) > 50000 else f"'{thumbnail[:50000]}'"},
        {'NULL' if not properties else f"'{properties}'"}, '{user_id}', {stability}, NOW());"""
        
        print(sql)
        print()


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Generate sample molecules for the library"
    )
    parser.add_argument(
        "--output",
        type=str,
        default="sample_molecules.json",
        help="Output JSON file path"
    )
    parser.add_argument(
        "--sql",
        action="store_true",
        help="Print Supabase INSERT SQL statements"
    )
    parser.add_argument(
        "--user-id",
        type=str,
        default="YOUR_USER_ID",
        help="User ID for SQL generation"
    )
    
    args = parser.parse_args()
    
    # Generate molecules
    molecules = save_to_json(args.output)
    
    # Print SQL if requested
    if args.sql:
        print_supabase_insert_sql(molecules, args.user_id)
    
    print(f"\n📊 Summary:")
    print(f"  • Total molecules: {len(molecules)}")
    print(f"  • With molfiles: {sum(1 for m in molecules if m['molfile'])}")
    print(f"  • With thumbnails: {sum(1 for m in molecules if m['thumbnail_b64'])}")
    print(f"\n💡 Next steps:")
    print(f"  1. Review {args.output}")
    print(f"  2. Use the JSON to seed your database")
    print(f"  3. Or use the SQL statements (with --sql flag)")


if __name__ == "__main__":
    main()

