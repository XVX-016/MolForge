import type { Atom, Bond, Molecule } from '../../types/molecule'
import { VALENCE_LIMITS } from '../bondRules'

export interface ValidationIssue {
  type: 'overbonded' | 'underbonded' | 'invalid_element'
  atomId: string
  element: string
  currentBonds: number
  maxBonds: number
  message: string
}

/**
 * Validate molecule structure and return issues
 */
export function validateStructure(mol: Molecule): ValidationIssue[] {
  const issues: ValidationIssue[] = []

  mol.atoms.forEach(atom => {
    const bondCount = mol.bonds.filter(
      b => b.atom1 === atom.id || b.atom2 === atom.id
    ).length

    const maxBonds = VALENCE_LIMITS[atom.element] || 4

    if (bondCount > maxBonds) {
      issues.push({
        type: 'overbonded',
        atomId: atom.id,
        element: atom.element,
        currentBonds: bondCount,
        maxBonds,
        message: `${atom.element} atom has ${bondCount} bonds, exceeds limit of ${maxBonds}`,
      })
    } else if (bondCount === 0 && atom.element !== 'H') {
      // Lone atoms (except H) might be underbonded
      issues.push({
        type: 'underbonded',
        atomId: atom.id,
        element: atom.element,
        currentBonds: bondCount,
        maxBonds,
        message: `${atom.element} atom has no bonds`,
      })
    }
  })

  return issues
}

