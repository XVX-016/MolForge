/**
 * SMILES (Simplified Molecular Input Line Entry System) utilities
 * 
 * SMILES is a string format that describes the structure of a molecule.
 * 
 * Examples:
 * - Water: O
 * - Ethanol: CCO
 * - Benzene: c1ccccc1
 * - Aspirin: CC(=O)Oc1ccccc1C(=O)O
 */

/**
 * Validate if a string is a valid SMILES format
 * This is a basic check - for production, use a proper SMILES parser
 */
export function isValidSMILES(smiles: string): boolean {
  if (!smiles || typeof smiles !== 'string') return false;
  
  // Basic validation: SMILES should contain only valid characters
  // This is a simplified check - real validation requires parsing
  const validChars = /^[A-Za-z0-9\[\]()=#@+\-.,;:*/\\$%]+$/;
  return validChars.test(smiles.trim()) && smiles.trim().length > 0;
}

/**
 * Normalize SMILES string (remove whitespace, convert to canonical form if possible)
 */
export function normalizeSMILES(smiles: string): string {
  return smiles.trim();
}

/**
 * Extract molecular formula from SMILES (basic implementation)
 * For production, use a proper SMILES parser like RDKit
 */
export function smilesToFormula(smiles: string): string | null {
  // This is a placeholder - real implementation requires SMILES parsing
  // For now, return null to indicate formula needs to be calculated elsewhere
  return null;
}

/**
 * Check if two SMILES strings represent the same molecule
 * Note: This requires canonicalization which is complex
 * For production, use a library like RDKit or OpenEye
 */
export function areSMILESEqual(smiles1: string, smiles2: string): boolean {
  // Basic comparison - for production, use canonical SMILES
  return normalizeSMILES(smiles1) === normalizeSMILES(smiles2);
}

/**
 * Search molecules by SMILES pattern (basic substring search)
 * For advanced search, use RDKit substructure matching
 */
export function searchBySMILES(molecules: Array<{ smiles?: string }>, query: string): Array<{ smiles?: string }> {
  const normalizedQuery = normalizeSMILES(query).toLowerCase();
  return molecules.filter((mol) => {
    if (!mol.smiles) return false;
    return normalizeSMILES(mol.smiles).toLowerCase().includes(normalizedQuery);
  });
}

