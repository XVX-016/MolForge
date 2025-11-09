import { describe, it, expect, beforeEach } from 'vitest'
import { MoleculeGraph } from '@biosynth/engine'
import { useMoleculeStore } from '../moleculeStore'
import { useHistoryStore, pushState, undo, redo, clearHistory } from '../historyStore'

describe('HistoryStore', () => {
  beforeEach(() => {
    useMoleculeStore.getState().reset()
    clearHistory()
  })

  it('pushes state to history', () => {
    const molecule = new MoleculeGraph()
    molecule.addAtom({ element: 'C', position: [0, 0, 0] })
    
    useMoleculeStore.getState().setMolecule(molecule)
    pushState()
    
    expect(useHistoryStore.getState().canUndo).toBe(false) // First state, can't undo
  })

  it('undoes and redoes actions', () => {
    const molecule1 = new MoleculeGraph()
    molecule1.addAtom({ element: 'C', position: [0, 0, 0] })
    useMoleculeStore.getState().setMolecule(molecule1)
    pushState()
    
    const molecule2 = molecule1.clone()
    molecule2.addAtom({ element: 'H', position: [1, 0, 0] })
    useMoleculeStore.getState().setMolecule(molecule2)
    pushState()
    
    expect(useMoleculeStore.getState().currentMolecule?.atoms.size).toBe(2)
    
    undo()
    expect(useMoleculeStore.getState().currentMolecule?.atoms.size).toBe(1)
    
    redo()
    expect(useMoleculeStore.getState().currentMolecule?.atoms.size).toBe(2)
  })

  it('tracks canUndo and canRedo', () => {
    const molecule = new MoleculeGraph()
    molecule.addAtom({ element: 'C', position: [0, 0, 0] })
    
    useMoleculeStore.getState().setMolecule(molecule)
    pushState()
    
    const molecule2 = molecule.clone()
    molecule2.addAtom({ element: 'H', position: [1, 0, 0] })
    useMoleculeStore.getState().setMolecule(molecule2)
    pushState()
    
    expect(useHistoryStore.getState().canUndo).toBe(true)
    expect(useHistoryStore.getState().canRedo).toBe(false)
    
    undo()
    expect(useHistoryStore.getState().canUndo).toBe(false)
    expect(useHistoryStore.getState().canRedo).toBe(true)
  })
})

