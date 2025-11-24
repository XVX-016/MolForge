/**
 * Bond distance thresholds and standard bond lengths
 * Used for auto-bonding logic
 */

export interface BondRule {
  threshold: number // Maximum distance for auto-bonding (Å)
  standardLength: number // Standard bond length (Å)
}

export const BOND_RULES: Record<string, BondRule> = {
  'H-H': { threshold: 0.8, standardLength: 0.74 },
  'C-C': { threshold: 1.8, standardLength: 1.54 },
  'C-H': { threshold: 1.3, standardLength: 1.09 },
  'C-O': { threshold: 1.6, standardLength: 1.43 },
  'C-N': { threshold: 1.7, standardLength: 1.47 },
  'O-H': { threshold: 1.2, standardLength: 0.96 },
  'N-H': { threshold: 1.3, standardLength: 1.01 },
  'O-O': { threshold: 1.6, standardLength: 1.48 },
  'N-N': { threshold: 1.7, standardLength: 1.45 },
  'C-F': { threshold: 1.6, standardLength: 1.35 },
  'C-Cl': { threshold: 2.0, standardLength: 1.76 },
  'C-S': { threshold: 2.0, standardLength: 1.82 },
  'C-P': { threshold: 2.1, standardLength: 1.87 },
}

/**
 * Get bond rule for a pair of elements
 */
export function getBondRule(element1: string, element2: string): BondRule {
  const key1 = `${element1}-${element2}`
  const key2 = `${element2}-${element1}`
  return BOND_RULES[key1] || BOND_RULES[key2] || { threshold: 1.8, standardLength: 1.5 }
}

/**
 * Typical valence limits for elements
 */
export const VALENCE_LIMITS: Record<string, number> = {
  H: 1,
  C: 4,
  N: 3,
  O: 2,
  F: 1,
  Cl: 1,
  Br: 1,
  I: 1,
  S: 6,
  P: 5,
}

