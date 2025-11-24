import type { Tool } from './toolInterface'
import * as THREE from 'three'

let pendingAtomId: string | null = null

const bondTool: Tool = {
  name: 'bond',
  onPointerDown: (ev: any, store: any) => {
    const pickedAtomId = ev.object?.userData?.atomId
    if (!pickedAtomId) {
      // Clicked empty space - cancel bonding
      pendingAtomId = null
      store.setSelectedAtomId(null)
      return
    }
    
    if (!pendingAtomId) {
      // First atom selected
      pendingAtomId = pickedAtomId
      store.setSelectedAtomId(pickedAtomId)
    } else {
      // Second atom selected - create bond with proper length
      if (pendingAtomId !== pickedAtomId) {
        const atom1 = store.molecule.atoms.find((a: any) => a.id === pendingAtomId)
        const atom2 = store.molecule.atoms.find((a: any) => a.id === pickedAtomId)
        
        if (atom1 && atom2) {
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
          store.addBond(pendingAtomId, pickedAtomId, 1)
        }
      }
      
      pendingAtomId = null
      store.setSelectedAtomId(null)
    }
  },
  onPointerUp: () => {
    // Reset on pointer up if needed
  }
}

export default bondTool
