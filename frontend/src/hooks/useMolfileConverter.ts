import { useState, useCallback } from 'react';
import { convertSMILESToMolfile } from '../lib/api';
import { updateMolecule, type SupabaseMolecule } from '../lib/supabaseMoleculeStore';

/**
 * Hook to auto-convert SMILES to molfile when missing
 * Returns a function that can be called to convert and update a molecule
 */
export function useMolfileConverter(userId: string | null) {
  const [converting, setConverting] = useState<string | null>(null); // molecule ID being converted

  const convertAndUpdate = useCallback(
    async (molecule: SupabaseMolecule): Promise<string | null> => {
      // Skip if no SMILES or already has molfile
      if (!molecule.smiles || molecule.molfile || !molecule.id || !userId) {
        return molecule.molfile || null;
      }

      // Skip if already converting this molecule
      if (converting === molecule.id) {
        return null;
      }

      setConverting(molecule.id);

      try {
        // Convert SMILES to molfile
        const result = await convertSMILESToMolfile(molecule.smiles);
        
        // Update molecule in Supabase
        await updateMolecule(userId, molecule.id, {
          molfile: result.molfile,
        });

        return result.molfile;
      } catch (error) {
        console.error('Failed to convert SMILES to molfile:', error);
        return null;
      } finally {
        setConverting(null);
      }
    },
    [userId, converting]
  );

  return { convertAndUpdate, converting };
}

