/**
 * 3D Service
 * 
 * Phase 9: 3D Preview Integration
 * 
 * Provides 3D coordinate generation using RDKit backend:
 * - ETKDG embedding
 * - Geometry optimization
 * - Fallback handling
 */

import type { Molecule } from '../Molecule'

export interface ThreeDOptions {
  method?: 'etkdg' | 'embed'
  optimize?: boolean
}

/**
 * Generate 3D coordinates for molecule
 */
export async function generate3DCoordinates(
  molecule: Molecule,
  options: ThreeDOptions = {}
): Promise<Molecule> {
  if (molecule.isEmpty()) {
    return molecule
  }

  // Convert molecule to backend format
  const moleculeData = moleculeToBackendFormat(molecule)

  try {
    const response = await fetch('/api/molecule/generate-3d', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        molecule: moleculeData,
        method: options.method || 'etkdg',
        optimize: options.optimize !== false, // Default to true
      }),
    })

    if (!response.ok) {
      throw new Error(`Failed to generate 3D coordinates: ${response.statusText}`)
    }

    const result = await response.json()
    return moleculeFromBackendFormat(result.molecule)
  } catch (error) {
    console.error('Error generating 3D coordinates:', error)
    // Fallback: return original molecule
    return molecule
  }
}

/**
 * Generate 3D coordinates from SMILES
 */
export async function generate3DFromSMILES(
  smiles: string,
  options: ThreeDOptions = {}
): Promise<Molecule> {
  try {
    const response = await fetch('/api/molecule/generate-3d', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        smiles,
        method: options.method || 'etkdg',
        optimize: options.optimize !== false,
      }),
    })

    if (!response.ok) {
      throw new Error(`Failed to generate 3D coordinates: ${response.statusText}`)
    }

    const result = await response.json()
    return moleculeFromBackendFormat(result.molecule)
  } catch (error) {
    console.error('Error generating 3D coordinates from SMILES:', error)
    throw error
  }
}

/**
 * Convert molecule to backend format
 */
function moleculeToBackendFormat(molecule: Molecule): any {
  const atoms = molecule.getAtoms().map(atom => ({
    id: atom.id,
    element: atom.element,
    position: atom.position,
    charge: atom.charge,
    formalCharge: atom.formalCharge,
  }))

  const bonds = molecule.getBonds().map(bond => ({
    id: bond.id,
    atom1: bond.atom1,
    atom2: bond.atom2,
    order: bond.order,
    type: bond.type,
  }))

  return {
    atoms,
    bonds,
    metadata: molecule.getMetadata(),
  }
}

/**
 * Convert backend format to Molecule
 */
function moleculeFromBackendFormat(data: any): Molecule {
  const { Molecule } = require('../Molecule')
  const molecule = new Molecule()

  // Add atoms
  data.atoms?.forEach((atom: any) => {
    molecule.addAtom({
      id: atom.id,
      element: atom.element,
      position: atom.position || [0, 0, 0],
      charge: atom.charge,
      formalCharge: atom.formalCharge,
    })
  })

  // Add bonds
  data.bonds?.forEach((bond: any) => {
    molecule.addBond({
      id: bond.id,
      atom1: bond.atom1,
      atom2: bond.atom2,
      order: bond.order || 1,
      type: bond.type,
    })
  })

  if (data.metadata) {
    molecule.setMetadata(data.metadata)
  }

  return molecule
}

