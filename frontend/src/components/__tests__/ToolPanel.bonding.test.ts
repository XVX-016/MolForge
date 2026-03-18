import { MoleculeGraph } from '@biosynth/engine'
import { beforeEach, describe, expect, it } from 'vitest'
import { createBondSafe } from '../../kernel/bonds'
import { getValence, shouldFormBond } from '../../kernel/valence'

interface EngineAtomLike {
  id: string
  element: string
  position: [number, number, number]
}

describe('ToolPanel -> kernel bond integration', () => {
  let molecule: MoleculeGraph

  beforeEach(() => {
    molecule = new MoleculeGraph()
  })

  it('creates valid bond between H and O using kernel function', () => {
    const hId = molecule.addAtom({ element: 'H', position: [0, 0, 0] })
    const oId = molecule.addAtom({ element: 'O', position: [0.9, 0, 0] })

    const bondId = createBondSafe(molecule, hId, oId)
    expect(bondId).toBeDefined()
    expect(bondId).not.toBeNull()
  })

  it('validates valence rules before creating bond', () => {
    const hId = molecule.addAtom({ element: 'H', position: [0, 0, 0] })
    const oId = molecule.addAtom({ element: 'O', position: [0.9, 0, 0] })

    const bond1 = createBondSafe(molecule, hId, oId)
    expect(bond1).not.toBeNull()

    const otherId = molecule.addAtom({ element: 'H', position: [2, 0, 0] })
    const bond2 = createBondSafe(molecule, hId, otherId)
    expect(bond2).toBeNull()
  })

  it('uses valence helper functions correctly', () => {
    expect(getValence('H')).toBe(1)
    expect(getValence('C')).toBe(4)
    expect(getValence('O')).toBe(2)

    const h: EngineAtomLike = { id: 'h1', element: 'H', position: [0, 0, 0] }
    const o: EngineAtomLike = { id: 'o1', element: 'O', position: [0.9, 0, 0] }

    expect(shouldFormBond(h, o, [], 1.2)).toBe(true)
  })

  it('prevents duplicate bonds', () => {
    const a1 = molecule.addAtom({ element: 'C', position: [0, 0, 0] })
    const a2 = molecule.addAtom({ element: 'C', position: [1.4, 0, 0] })

    const bond1 = createBondSafe(molecule, a1, a2)
    expect(bond1).not.toBeNull()

    const bond2 = createBondSafe(molecule, a1, a2)
    expect(bond2).toBe(bond1)
  })
})
