/**
 * Ring Templates
 * 
 * Phase 11: Toolbar Rebuild
 * 
 * Pre-defined ring structures for quick molecule building
 */

import type { Molecule } from '../Molecule'
import { nanoid } from 'nanoid'

export interface RingTemplate {
  id: string
  name: string
  smiles: string
  description: string
}

export const RING_TEMPLATES: RingTemplate[] = [
  {
    id: 'benzene',
    name: 'Benzene',
    smiles: 'c1ccccc1',
    description: '6-membered aromatic ring',
  },
  {
    id: 'cyclopropane',
    name: 'Cyclopropane',
    smiles: 'C1CC1',
    description: '3-membered ring',
  },
  {
    id: 'cyclobutane',
    name: 'Cyclobutane',
    smiles: 'C1CCC1',
    description: '4-membered ring',
  },
  {
    id: 'cyclopentane',
    name: 'Cyclopentane',
    smiles: 'C1CCCC1',
    description: '5-membered ring',
  },
  {
    id: 'cyclohexane',
    name: 'Cyclohexane',
    smiles: 'C1CCCCC1',
    description: '6-membered ring',
  },
  {
    id: 'pyridine',
    name: 'Pyridine',
    smiles: 'c1ccncc1',
    description: '6-membered aromatic with N',
  },
  {
    id: 'pyrrole',
    name: 'Pyrrole',
    smiles: 'c1cc[nH]c1',
    description: '5-membered aromatic with N',
  },
  {
    id: 'furan',
    name: 'Furan',
    smiles: 'c1ccoc1',
    description: '5-membered aromatic with O',
  },
  {
    id: 'thiophene',
    name: 'Thiophene',
    smiles: 'c1ccsc1',
    description: '5-membered aromatic with S',
  },
]

/**
 * Load ring template into molecule
 * Uses 2D layout generation to position atoms
 */
export async function loadRingTemplate(
  template: RingTemplate,
  centerX: number = 0,
  centerY: number = 0
): Promise<Molecule> {
  const { generate2DLayoutFromSMILES } = await import('../layout')
  
  try {
    const molecule = await generate2DLayoutFromSMILES(template.smiles)
    
    // Center the molecule at the specified position
    const atoms = molecule.getAtoms()
    if (atoms.length === 0) {
      return molecule
    }
    
    // Calculate centroid
    let sumX = 0
    let sumY = 0
    atoms.forEach(atom => {
      sumX += atom.position[0]
      sumY += atom.position[1]
    })
    const centroidX = sumX / atoms.length
    const centroidY = sumY / atoms.length
    
    // Translate to center position
    const offsetX = centerX - centroidX
    const offsetY = centerY - centroidY
    
    // Create new molecule with translated positions
    const { Molecule } = await import('../Molecule')
    const centeredMolecule = new Molecule()
    
    atoms.forEach(atom => {
      centeredMolecule.addAtom({
        id: atom.id,
        element: atom.element,
        position: [
          atom.position[0] + offsetX,
          atom.position[1] + offsetY,
          atom.position[2],
        ],
        charge: atom.charge,
        formalCharge: atom.formalCharge,
      })
    })
    
    molecule.getBonds().forEach(bond => {
      centeredMolecule.addBond({
        id: bond.id,
        atom1: bond.atom1,
        atom2: bond.atom2,
        order: bond.order,
        type: bond.type,
      })
    })
    
    return centeredMolecule
  } catch (error) {
    console.error(`Failed to load ring template ${template.name}:`, error)
    throw error
  }
}

