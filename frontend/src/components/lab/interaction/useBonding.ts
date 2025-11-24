import { useRef } from 'react'
import { useLabStore } from '../../../store/labStore'

/**
 * Hook for manual bond creation
 * Tracks pending atom selection for bond creation
 */
export function useBonding() {
  const pendingAtomRef = useRef<string | null>(null)
  const store = useLabStore.getState()

  const handleAtomClick = (atomId: string): boolean => {
    if (!pendingAtomRef.current) {
      // First atom selected
      pendingAtomRef.current = atomId
      store.setSelectedAtomId(atomId)
      return false // Bond not created yet
    } else {
      // Second atom selected - create bond
      if (pendingAtomRef.current !== atomId) {
        // Standard bond lengths (Å)
        const bondLengths: Record<string, number> = {
          'C-C': 1.54,
          'C-H': 1.09,
          'O-H': 0.96,
          'N-H': 1.01,
          'C-O': 1.43,
          'C-N': 1.47,
          'O-O': 1.48,
          'N-N': 1.45,
        }

        const atom1 = store.molecule.atoms.find(a => a.id === pendingAtomRef.current)
        const atom2 = store.molecule.atoms.find(a => a.id === atomId)
        
        if (atom1 && atom2) {
          // Get bond length
          const key1 = `${atom1.element}-${atom2.element}`
          const key2 = `${atom2.element}-${atom1.element}`
          const targetLength = bondLengths[key1] || bondLengths[key2] || 1.5

          // Adjust atom2 position to match bond length
          const dir = new THREE.Vector3(
            atom2.position[0] - atom1.position[0],
            atom2.position[1] - atom1.position[1],
            atom2.position[2] - atom1.position[2]
          )
          const currentLength = dir.length()
          if (currentLength > 0) {
            dir.normalize().multiplyScalar(targetLength)
            const newPos: [number, number, number] = [
              atom1.position[0] + dir.x,
              atom1.position[1] + dir.y,
              atom1.position[2] + dir.z,
            ]
            store.moveAtom(atom2.id, newPos)
          }

          // Create bond
          store.addBond(pendingAtomRef.current, atomId, 1)
        }
      }
      
      pendingAtomRef.current = null
      store.setSelectedAtomId(null)
      return true // Bond created
    }
  }

  const cancelBonding = () => {
    pendingAtomRef.current = null
    store.setSelectedAtomId(null)
  }

  return { handleAtomClick, cancelBonding, pendingAtomId: pendingAtomRef.current }
}

