/**
 * MoleculeEditor - Main molecule editor component
 * 
 * Phase 2: Centralized editor state and operations
 * 
 * This component manages:
 * - Molecule state (using Molecule class from Phase 1)
 * - Drawing operations (add atom, add bond, delete, move)
 * - Selection state
 * - Hover state
 * - Tool mode
 */

import React, { useState, useCallback, useRef } from 'react'
import { Molecule } from '@/lib/molecule'
import type { EditorTool, Atom, Bond } from '@/lib/molecule'
import { CanvasLayer } from './CanvasLayer'
import { nanoid } from 'nanoid'

interface MoleculeEditorProps {
  width: number
  height: number
  initialMolecule?: Molecule
  tool?: EditorTool
  onMoleculeChange?: (molecule: Molecule) => void
  onAtomSelect?: (atomId: string | null) => void
  onBondSelect?: (bondId: string | null) => void
}

export function MoleculeEditor({
  width,
  height,
  initialMolecule,
  tool = 'select',
  onMoleculeChange,
  onAtomSelect,
  onBondSelect,
}: MoleculeEditorProps) {
  const [molecule, setMolecule] = useState<Molecule>(initialMolecule || new Molecule())
  const [selectedAtomId, setSelectedAtomId] = useState<string | null>(null)
  const [selectedBondId, setSelectedBondId] = useState<string | null>(null)
  const [hoveredAtomId, setHoveredAtomId] = useState<string | null>(null)
  const [hoveredBondId, setHoveredBondId] = useState<string | null>(null)
  const [elementToAdd, setElementToAdd] = useState<string>('C')
  const [bondOrder, setBondOrder] = useState<number>(1)
  const [scale, setScale] = useState(1)
  const [offsetX, setOffsetX] = useState(0)
  const [offsetY, setOffsetY] = useState(0)

  // Update molecule and notify parent
  const updateMolecule = useCallback((newMolecule: Molecule) => {
    setMolecule(newMolecule)
    onMoleculeChange?.(newMolecule)
  }, [onMoleculeChange])

  // Add atom at position
  const addAtom = useCallback((x: number, y: number, element: string = elementToAdd) => {
    const newMolecule = molecule.clone()
    
    // Convert screen coordinates to world coordinates
    const worldX = (x - width / 2 - offsetX) / scale
    const worldY = (y - height / 2 - offsetY) / scale
    const worldZ = 0

    const atom: Atom = {
      id: nanoid(),
      element,
      position: [worldX, worldY, worldZ],
    }

    newMolecule.addAtom(atom)
    updateMolecule(newMolecule)
    return atom.id
  }, [molecule, width, height, scale, offsetX, offsetY, elementToAdd, updateMolecule])

  // Add bond between two atoms
  const addBond = useCallback((atom1Id: string, atom2Id: string, order: number = bondOrder) => {
    if (atom1Id === atom2Id) return

    const newMolecule = molecule.clone()
    
    // Check if bond already exists
    const existingBonds = newMolecule.getBonds()
    for (const bond of existingBonds) {
      if (
        (bond.atom1 === atom1Id && bond.atom2 === atom2Id) ||
        (bond.atom1 === atom2Id && bond.atom2 === atom1Id)
      ) {
        return // Bond already exists
      }
    }

    const bond: Bond = {
      id: nanoid(),
      atom1: atom1Id,
      atom2: atom2Id,
      order,
    }

    newMolecule.addBond(bond)
    updateMolecule(newMolecule)
    return bond.id
  }, [molecule, bondOrder, updateMolecule])

  // Delete atom (and all its bonds)
  const deleteAtom = useCallback((atomId: string) => {
    const newMolecule = molecule.clone()
    newMolecule.removeAtom(atomId)
    updateMolecule(newMolecule)
    
    if (selectedAtomId === atomId) {
      setSelectedAtomId(null)
      onAtomSelect?.(null)
    }
  }, [molecule, selectedAtomId, updateMolecule, onAtomSelect])

  // Delete bond
  const deleteBond = useCallback((bondId: string) => {
    const newMolecule = molecule.clone()
    newMolecule.removeBond(bondId)
    updateMolecule(newMolecule)
    
    if (selectedBondId === bondId) {
      setSelectedBondId(null)
      onBondSelect?.(null)
    }
  }, [molecule, selectedBondId, updateMolecule, onBondSelect])

  // Move atom
  const moveAtom = useCallback((atomId: string, x: number, y: number) => {
    const newMolecule = molecule.clone()
    const worldX = (x - width / 2 - offsetX) / scale
    const worldY = (y - height / 2 - offsetY) / scale
    const atom = newMolecule.getAtom(atomId)
    if (atom) {
      newMolecule.updateAtomPosition(atomId, [worldX, worldY, atom.position[2]])
      updateMolecule(newMolecule)
    }
  }, [molecule, width, height, scale, offsetX, offsetY, updateMolecule])

  // Handle atom click
  const handleAtomClick = useCallback((atomId: string, event: MouseEvent) => {
    if (tool === 'select') {
      setSelectedAtomId(atomId)
      setSelectedBondId(null)
      onAtomSelect?.(atomId)
      onBondSelect?.(null)
    } else if (tool === 'bond' && selectedAtomId && selectedAtomId !== atomId) {
      // Create bond between selected atom and clicked atom
      addBond(selectedAtomId, atomId)
      setSelectedAtomId(null)
      onAtomSelect?.(null)
    } else if (tool === 'delete') {
      deleteAtom(atomId)
    } else if (tool === 'add-atom') {
      // Place atom at click position (handled by canvas click)
    }
  }, [tool, selectedAtomId, addBond, deleteAtom, onAtomSelect, onBondSelect])

  // Handle bond click
  const handleBondClick = useCallback((bondId: string, event: MouseEvent) => {
    if (tool === 'select') {
      setSelectedBondId(bondId)
      setSelectedAtomId(null)
      onBondSelect?.(bondId)
      onAtomSelect?.(null)
    } else if (tool === 'delete') {
      deleteBond(bondId)
    }
  }, [tool, deleteBond, onBondSelect, onAtomSelect])

  // Handle canvas click (for adding atoms)
  const handleCanvasClick = useCallback((x: number, y: number) => {
    if (tool === 'add-atom') {
      addAtom(x, y)
    }
  }, [tool, addAtom])

  // Clear molecule
  const clear = useCallback(() => {
    const newMolecule = new Molecule()
    updateMolecule(newMolecule)
    setSelectedAtomId(null)
    setSelectedBondId(null)
    onAtomSelect?.(null)
    onBondSelect?.(null)
  }, [updateMolecule, onAtomSelect, onBondSelect])

  return (
    <div className="molecule-editor" style={{ width, height, position: 'relative' }}>
      <CanvasLayer
        molecule={molecule}
        selectedAtomId={selectedAtomId}
        selectedBondId={selectedBondId}
        hoveredAtomId={hoveredAtomId}
        hoveredBondId={hoveredBondId}
        width={width}
        height={height}
        scale={scale}
        offsetX={offsetX}
        offsetY={offsetY}
        onAtomClick={handleAtomClick}
        onBondClick={handleBondClick}
        onAtomHover={setHoveredAtomId}
        onBondHover={setHoveredBondId}
      />
    </div>
  )
}

// Export editor operations for external use
export type MoleculeEditorRef = {
  addAtom: (x: number, y: number, element?: string) => string
  addBond: (atom1Id: string, atom2Id: string, order?: number) => string | null
  deleteAtom: (atomId: string) => void
  deleteBond: (bondId: string) => void
  moveAtom: (atomId: string, x: number, y: number) => void
  clear: () => void
  getMolecule: () => Molecule
}

