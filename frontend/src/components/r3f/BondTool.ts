import { useEffect, useRef } from 'react'
import { useMoleculeStore } from '../../store/moleculeStore'
import { addBond } from '../../lib/engineAdapter'

/**
 * BondTool - Handles bond creation between selected atoms
 * Creates a bond when two atoms are selected within 2 seconds
 */
export function useBondTool() {
  const selectedAtomId = useMoleculeStore((state) => state.selectedAtomId)
  const currentMolecule = useMoleculeStore((state) => state.currentMolecule)
  const tool = useMoleculeStore((state) => state.tool)
  const currentBondOrder = useMoleculeStore((state) => state.currentBondOrder)
  const firstSelectedRef = useRef<string | null>(null)
  const timeoutRef = useRef<NodeJS.Timeout | null>(null)

  useEffect(() => {
    // Only work when bond tool is active
    if (tool !== 'bond') {
      firstSelectedRef.current = null
      if (timeoutRef.current) {
        clearTimeout(timeoutRef.current)
        timeoutRef.current = null
      }
      return
    }

    // Reset timeout when selection changes
    if (timeoutRef.current) {
      clearTimeout(timeoutRef.current)
      timeoutRef.current = null
    }

    if (!selectedAtomId || !currentMolecule) {
      firstSelectedRef.current = null
      return
    }

    // If no first atom selected, store it
    if (!firstSelectedRef.current) {
      firstSelectedRef.current = selectedAtomId
      return
    }

    // If same atom selected, ignore
    if (firstSelectedRef.current === selectedAtomId) {
      return
    }

    // Two different atoms selected - create bond
    const atom1 = firstSelectedRef.current
    const atom2 = selectedAtomId

    // Check if bond already exists
    const existingBond = Array.from(currentMolecule.bonds.values()).find(
      (b) =>
        (b.a1 === atom1 && b.a2 === atom2) ||
        (b.a1 === atom2 && b.a2 === atom1)
    )

    if (!existingBond) {
      // Create bond with current bond order
      addBond(atom1, atom2, currentBondOrder)
    }

    // Reset after bond creation
    firstSelectedRef.current = null
  }, [selectedAtomId, currentMolecule, tool, currentBondOrder])

  // Handle Escape key to cancel bond creation
  useEffect(() => {
    const handleEscape = (e: KeyboardEvent) => {
      if (e.key === 'Escape' && tool === 'bond') {
        firstSelectedRef.current = null
        if (timeoutRef.current) {
          clearTimeout(timeoutRef.current)
          timeoutRef.current = null
        }
      }
    }

    window.addEventListener('keydown', handleEscape)
    return () => {
      window.removeEventListener('keydown', handleEscape)
    }
  }, [tool])

  return {
    firstSelectedAtomId: firstSelectedRef.current,
    isWaitingForSecondAtom: firstSelectedRef.current !== null,
  }
}

