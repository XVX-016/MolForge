#!/usr/bin/env python3
"""
Batch molecule asset generation script
Generates 2D thumbnails and 3D molfiles from molecule files

Usage:
    python generate_molecule_assets.py --input-dir molecules_raw --output-dir public/molecules
"""
import os
import argparse
import json
from pathlib import Path
from typing import Optional
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from backend.services.thumbnail_service import ThumbnailService


def load_molecule_from_file(file_path: Path) -> Optional[Chem.Mol]:
    """Load molecule from various file formats"""
    suffix = file_path.suffix.lower()
    
    try:
        if suffix == '.mol':
            return Chem.MolFromMolFile(str(file_path))
        elif suffix == '.sdf':
            # SDF files can contain multiple molecules, take first
            supplier = Chem.SDMolSupplier(str(file_path))
            return next(supplier, None)
        elif suffix == '.pdb':
            return Chem.MolFromPDBFile(str(file_path))
        elif suffix == '.xyz':
            # XYZ format support (basic)
            with open(file_path, 'r') as f:
                lines = f.readlines()
            # Simple XYZ parser - you may want to use a proper library
            return None  # Placeholder - implement XYZ parsing if needed
        elif suffix == '.smi':
            # SMILES file
            with open(file_path, 'r') as f:
                smiles = f.readline().strip()
            return Chem.MolFromSmiles(smiles)
        else:
            print(f"⚠️  Unsupported file format: {suffix}")
            return None
    except Exception as e:
        print(f"❌ Error reading {file_path}: {e}")
        return None


def process_molecule_file(
    input_file: Path,
    output_dir: Path,
    generate_3d: bool = True
) -> dict:
    """
    Process a single molecule file and generate assets
    
    Returns:
        Dictionary with generated file paths and metadata
    """
    base_name = input_file.stem
    
    # Load molecule
    mol = load_molecule_from_file(input_file)
    if mol is None:
        return {"error": f"Could not load molecule from {input_file.name}"}
    
    results = {
        "name": base_name,
        "input_file": str(input_file),
        "generated_files": []
    }
    
    try:
        # Generate 2D thumbnail
        thumbnail_path = output_dir / f"{base_name}.png"
        img_bytes = ThumbnailService.generate_2d_thumbnail(
            molfile=Chem.MolToMolBlock(mol),
            size=(600, 600)
        )
        thumbnail_path.write_bytes(img_bytes)
        results["generated_files"].append(str(thumbnail_path))
        print(f"  ✓ Generated 2D thumbnail: {thumbnail_path.name}")
        
        # Generate 3D structure if requested
        if generate_3d:
            try:
                mol_3d = ThumbnailService.generate_3d_structure(
                    molfile=Chem.MolToMolBlock(mol)
                )
                molfile_3d = Chem.MolToMolBlock(mol_3d)
                
                # Save 3D molfile
                molfile_path = output_dir / f"{base_name}_3d.mol"
                molfile_path.write_text(molfile_3d)
                results["generated_files"].append(str(molfile_path))
                print(f"  ✓ Generated 3D molfile: {molfile_path.name}")
                
                # Extract SMILES if available
                smiles = Chem.MolToSmiles(mol)
                if smiles:
                    results["smiles"] = smiles
                    
            except Exception as e:
                print(f"  ⚠️  Could not generate 3D structure: {e}")
                results["3d_error"] = str(e)
        
        # Extract metadata
        try:
            smiles = Chem.MolToSmiles(mol)
            if smiles:
                results["smiles"] = smiles
        except Exception:
            pass
        
        # Save metadata JSON
        metadata_path = output_dir / f"{base_name}_metadata.json"
        metadata_path.write_text(json.dumps(results, indent=2))
        results["generated_files"].append(str(metadata_path))
        
    except Exception as e:
        results["error"] = str(e)
        print(f"  ❌ Error processing {input_file.name}: {e}")
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Generate molecule assets (thumbnails, 3D molfiles) from raw molecule files"
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        default="molecules_raw",
        help="Input directory containing molecule files (.mol, .sdf, .pdb, etc.)"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="public/molecules",
        help="Output directory for generated assets"
    )
    parser.add_argument(
        "--no-3d",
        action="store_true",
        help="Skip 3D structure generation (faster)"
    )
    parser.add_argument(
        "--extensions",
        nargs="+",
        default=[".mol", ".sdf", ".pdb", ".smi"],
        help="File extensions to process"
    )
    
    args = parser.parse_args()
    
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if not input_dir.exists():
        print(f"❌ Input directory does not exist: {input_dir}")
        return
    
    # Find all molecule files
    molecule_files = []
    for ext in args.extensions:
        molecule_files.extend(input_dir.glob(f"*{ext}"))
        molecule_files.extend(input_dir.glob(f"*{ext.upper()}"))
    
    if not molecule_files:
        print(f"⚠️  No molecule files found in {input_dir}")
        print(f"   Looking for extensions: {', '.join(args.extensions)}")
        return
    
    print(f"📦 Found {len(molecule_files)} molecule file(s)")
    print(f"📁 Output directory: {output_dir}")
    print(f"{'=' * 60}\n")
    
    # Process each file
    results = []
    for i, mol_file in enumerate(molecule_files, 1):
        print(f"[{i}/{len(molecule_files)}] Processing {mol_file.name}...")
        result = process_molecule_file(
            mol_file,
            output_dir,
            generate_3d=not args.no_3d
        )
        results.append(result)
        print()
    
    # Summary
    successful = sum(1 for r in results if "error" not in r)
    print(f"{'=' * 60}")
    print(f"✅ Successfully processed: {successful}/{len(results)}")
    print(f"📁 Output directory: {output_dir}")
    print(f"\nGenerated files:")
    for result in results:
        if "generated_files" in result:
            for file_path in result["generated_files"]:
                print(f"  • {Path(file_path).name}")


if __name__ == "__main__":
    main()

