import { describe, it, expect, beforeEach } from 'vitest'
import { MoleculeGraph } from '@biosynth/engine'
import { useMoleculeStore } from '../../store/moleculeStore'
import { addAtom, removeAtom, removeBond } from '../engineAdapter'

describe('engineAdapter - Atom Operations', () => {
  beforeEach(() => {
    useMoleculeStore.getState().reset()
  })

  it('adds atom to molecule', () => {
    addAtom('C', [0, 0, 0])
    
    const molecule = useMoleculeStore.getState().currentMolecule
    expect(molecule).not.toBeNull()
    expect(molecule?.atoms.size).toBe(1)
    
    const atom = Array.from(molecule!.atoms.values())[0]
    expect(atom.element).toBe('C')
    expect(atom.position).toEqual([0, 0, 0])
  })

  it('adds atom to existing molecule', () => {
    const molecule = new MoleculeGraph()
    molecule.addAtom({ element: 'C', position: [0, 0, 0] })
    useMoleculeStore.getState().setMolecule(molecule)
    
    addAtom('H', [1, 0, 0])
    
    const updated = useMoleculeStore.getState().currentMolecule
    expect(updated?.atoms.size).toBe(2)
  })

  it('removes atom from molecule', () => {
    const molecule = new MoleculeGraph()
    const atomId = molecule.addAtom({ element: 'C', position: [0, 0, 0] })
    useMoleculeStore.getState().setMolecule(molecule)
    
    removeAtom(atomId)
    
    const updated = useMoleculeStore.getState().currentMolecule
    expect(updated?.atoms.size).toBe(0)
  })

  it('removes bond from molecule', () => {
    const molecule = new MoleculeGraph()
    const atom1Id = molecule.addAtom({ element: 'C', position: [0, 0, 0] })
    const atom2Id = molecule.addAtom({ element: 'H', position: [1, 0, 0] })
    const bondId = molecule.addBond(atom1Id, atom2Id, 1)
    useMoleculeStore.getState().setMolecule(molecule)
    
    if (bondId) {
      removeBond(bondId)
      
      const updated = useMoleculeStore.getState().currentMolecule
      expect(updated?.bonds.size).toBe(0)
    }
  })
})

