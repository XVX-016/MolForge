"""
Thumbnail Generation Service
Generates 2D molecule thumbnails using RDKit
"""
import base64
import io
from typing import Optional, Tuple
from rdkit import Chem
from rdkit.Chem import Draw, AllChem


class ThumbnailService:
    """Service for generating molecule thumbnails"""
    
    @staticmethod
    def generate_2d_thumbnail(
        molfile: Optional[str] = None,
        smiles: Optional[str] = None,
        size: Tuple[int, int] = (600, 600),
        format: str = 'PNG'
    ) -> bytes:
        """
        Generate a 2D thumbnail image from molfile or SMILES
        
        Args:
            molfile: V2000/V3000 molfile string
            smiles: SMILES string
            size: Image size (width, height)
            format: Image format ('PNG' or 'SVG')
            
        Returns:
            Image bytes
        """
        # Load molecule from molfile or SMILES
        mol = None
        if molfile:
            mol = Chem.MolFromMolBlock(molfile)
        elif smiles:
            mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            raise ValueError("Could not parse molecule from provided molfile or SMILES")
        
        # Generate 2D coordinates if not present
        try:
            AllChem.Compute2DCoords(mol)
        except Exception:
            # If 2D coords already exist or computation fails, continue
            pass
        
        # Generate image
        if format.upper() == 'SVG':
            img_data = Draw.MolToSVG(mol, size=size)
            return img_data.encode('utf-8')
        else:
            img = Draw.MolToImage(mol, size=size)
            img_bytes = io.BytesIO()
            img.save(img_bytes, format='PNG')
            return img_bytes.getvalue()
    
    @staticmethod
    def generate_2d_thumbnail_base64(
        molfile: Optional[str] = None,
        smiles: Optional[str] = None,
        size: tuple[int, int] = (600, 600)
    ) -> str:
        """
        Generate a 2D thumbnail and return as base64 string
        
        Args:
            molfile: V2000/V3000 molfile string
            smiles: SMILES string
            size: Image size (width, height)
            
        Returns:
            Base64-encoded image string (data URI format)
        """
        img_bytes = ThumbnailService.generate_2d_thumbnail(molfile, smiles, size)
        img_b64 = base64.b64encode(img_bytes).decode('utf-8')
        return f"data:image/png;base64,{img_b64}"
    
    @staticmethod
    def generate_3d_structure(molfile: Optional[str] = None, smiles: Optional[str] = None):
        """
        Generate optimized 3D structure from molfile or SMILES
        
        Args:
            molfile: V2000/V3000 molfile string
            smiles: SMILES string
            
        Returns:
            RDKit Mol object with 3D coordinates
        """
        # Load molecule
        mol = None
        if molfile:
            mol = Chem.MolFromMolBlock(molfile)
        elif smiles:
            mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            raise ValueError("Could not parse molecule from provided molfile or SMILES")
        
        # Add hydrogens
        mol = Chem.AddHs(mol)
        
        # Generate 3D coordinates
        embed_result = AllChem.EmbedMolecule(mol, randomSeed=42)
        if embed_result == -1:
            # Try alternative method
            embed_result = AllChem.EmbedMolecule(mol, useRandomCoords=True)
            if embed_result == -1:
                raise ValueError("Failed to generate 3D coordinates for molecule")
        
        # Optimize geometry
        try:
            AllChem.MMFFOptimizeMolecule(mol)
        except Exception:
            try:
                AllChem.UFFOptimizeMolecule(mol)
            except Exception:
                # Continue with unoptimized coordinates
                pass
        
        return mol
    
    @staticmethod
    def generate_3d_molfile(
        molfile: Optional[str] = None,
        smiles: Optional[str] = None
    ) -> str:
        """
        Generate optimized 3D molfile from molfile or SMILES
        
        Args:
            molfile: V2000/V3000 molfile string
            smiles: SMILES string
            
        Returns:
            V2000 molfile string with 3D coordinates
        """
        mol = ThumbnailService.generate_3d_structure(molfile, smiles)
        return Chem.MolToMolBlock(mol)

