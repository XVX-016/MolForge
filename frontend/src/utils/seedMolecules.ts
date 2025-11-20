/**
 * Utility to seed the library with sample molecules
 * Can be used for testing or initial data population
 */

import { saveMolecule, type SupabaseMolecule } from '../lib/supabaseMoleculeStore';
import { convertSMILESToMolfile, generateThumbnailBase64 } from '../lib/api';

export interface SampleMolecule {
  name: string;
  smiles: string;
  formula: string;
  description?: string;
}

// Sample molecules for seeding
export const SAMPLE_MOLECULES: SampleMolecule[] = [
  {
    name: 'Aspirin',
    smiles: 'CC(=O)OC1=CC=CC=C1C(=O)O',
    formula: 'C9H8O4',
    description: 'Common pain reliever and anti-inflammatory drug',
  },
  {
    name: 'Caffeine',
    smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
    formula: 'C8H10N4O2',
    description: 'Stimulant found in coffee and tea',
  },
  {
    name: 'Glucose',
    smiles: 'C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O',
    formula: 'C6H12O6',
    description: 'Simple sugar, primary energy source',
  },
  {
    name: 'Ethanol',
    smiles: 'CCO',
    formula: 'C2H6O',
    description: 'Alcohol found in beverages',
  },
  {
    name: 'Paracetamol',
    smiles: 'CC(=O)NC1=CC=C(C=C1)O',
    formula: 'C8H9NO2',
    description: 'Pain reliever and fever reducer',
  },
  {
    name: 'Ibuprofen',
    smiles: 'CC(C)Cc1ccc(C(C)C(=O)O)cc1',
    formula: 'C13H18O2',
    description: 'Nonsteroidal anti-inflammatory drug',
  },
  {
    name: 'Vitamin C',
    smiles: 'C([C@@H]([C@H]1C(=C(C(=O)O1)O)O)O)O',
    formula: 'C6H8O6',
    description: 'Ascorbic acid, essential vitamin',
  },
  {
    name: 'Nicotine',
    smiles: 'CN1CCC[C@H]1c2cccnc2',
    formula: 'C10H14N2',
    description: 'Stimulant found in tobacco',
  },
  {
    name: 'Serotonin',
    smiles: 'C1=CC2=C(C=C1O)C(=CN2)CCN',
    formula: 'C10H12N2O',
    description: 'Neurotransmitter',
  },
  {
    name: 'Dopamine',
    smiles: 'C1=CC(=C(C=C1CCN)O)O',
    formula: 'C8H11NO2',
    description: 'Neurotransmitter',
  },
];

/**
 * Generate complete molecule data with molfile and thumbnail
 */
export async function generateMoleculeData(
  molecule: SampleMolecule
): Promise<Omit<SupabaseMolecule, 'id' | 'user_id' | 'created_at'>> {
  let molfile: string | undefined;
  let thumbnail_b64: string | undefined;

  try {
    // Generate 3D molfile
    const molfileResult = await convertSMILESToMolfile(molecule.smiles);
    molfile = molfileResult.molfile;
  } catch (error) {
    console.warn(`Failed to generate molfile for ${molecule.name}:`, error);
  }

  try {
    // Generate thumbnail
    const thumbResult = await generateThumbnailBase64({ smiles: molecule.smiles });
    thumbnail_b64 = thumbResult.thumbnail_b64;
  } catch (error) {
    console.warn(`Failed to generate thumbnail for ${molecule.name}:`, error);
  }

  // Generate mock properties
  const properties = {
    stability: 0.75 + (molecule.smiles.length % 25) / 100,
    toxicity: 0.1 + (molecule.smiles.length % 30) / 100,
    solubility: 0.5 + (molecule.smiles.length % 40) / 100,
    bioavailability: 0.6 + (molecule.smiles.length % 30) / 100,
    novelty: 0.3 + (molecule.smiles.length % 50) / 100,
  };

  return {
    name: molecule.name,
    smiles: molecule.smiles,
    formula: molecule.formula,
    molfile: molfile || undefined,
    thumbnail_b64: thumbnail_b64 || undefined,
    properties: JSON.stringify(properties),
  };
}

/**
 * Seed the library with sample molecules
 * @param userId - User ID to seed molecules for
 * @param onProgress - Optional callback for progress updates
 */
export async function seedMolecules(
  userId: string,
  onProgress?: (current: number, total: number) => void
): Promise<{ success: number; failed: number }> {
  if (!userId) {
    throw new Error('User ID is required to seed molecules');
  }

  let success = 0;
  let failed = 0;
  const total = SAMPLE_MOLECULES.length;

  console.log(`Seeding ${total} molecules for user ${userId}...`);

  for (let i = 0; i < SAMPLE_MOLECULES.length; i++) {
    const molecule = SAMPLE_MOLECULES[i];
    
    // Update progress
    if (onProgress) {
      onProgress(i + 1, total);
    }

    try {
      const moleculeData = await generateMoleculeData(molecule);
      await saveMolecule(userId, moleculeData);
      console.log(`✓ Saved: ${molecule.name}`);
      success++;
      
      // Small delay to avoid overwhelming the API
      await new Promise((resolve) => setTimeout(resolve, 200));
    } catch (error) {
      console.error(`✗ Failed to save ${molecule.name}:`, error);
      failed++;
    }
  }

  console.log(`\n✅ Seeding complete: ${success} succeeded, ${failed} failed`);
  return { success, failed };
}

