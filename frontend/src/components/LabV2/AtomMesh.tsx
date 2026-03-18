import { useCallback, useRef, useState } from 'react'
import type { Mesh } from 'three'
import type { ThreeEvent } from '@react-three/fiber'
import { useLabStore } from '../../store/labStore'
import type { Atom } from '../../types/molecule'
import { getAtomColor } from '../../utils/atomColors'
import { getElementSpec } from '../../utils/elements'

export default function AtomMesh({ atom }: { atom: Atom }) {
  const mesh = useRef<Mesh | null>(null)
  const { currentTool, addBond, deleteAtom, setSelectedAtomId, selectedAtomId, bondOrder } = useLabStore()
  const [hovered, setHovered] = useState(false)

  const color = getAtomColor(atom.element)
  const spec = getElementSpec(atom.element)
  const radius = spec ? spec.radius * 0.25 : atom.element === 'H' ? 0.18 : 0.28

  const handleClick = useCallback(
    (event: ThreeEvent<MouseEvent>) => {
      event.stopPropagation()

      if (currentTool === 'bond') {
        if (selectedAtomId && selectedAtomId !== atom.id) {
          addBond(selectedAtomId, atom.id, bondOrder)
          setSelectedAtomId(null)
        } else {
          setSelectedAtomId(atom.id)
        }
      } else if (currentTool === 'delete') {
        deleteAtom(atom.id)
      } else {
        setSelectedAtomId(atom.id)
      }
    },
    [addBond, atom.id, bondOrder, currentTool, deleteAtom, selectedAtomId, setSelectedAtomId]
  )

  const { x, y, z } = atom.position

  return (
    <mesh
      ref={mesh}
      position={[x, y, z]}
      onPointerOver={(event) => {
        event.stopPropagation()
        setHovered(true)
      }}
      onPointerOut={(event) => {
        event.stopPropagation()
        setHovered(false)
      }}
      onClick={handleClick}
      scale={hovered || selectedAtomId === atom.id ? 1.2 : 1}
    >
      <sphereGeometry args={[radius, 32, 32]} />
      <meshStandardMaterial color={color} metalness={0.1} roughness={0.4} />
      {(hovered || selectedAtomId === atom.id) && (
        <mesh>
          <sphereGeometry args={[radius * 1.22, 16, 16]} />
          <meshBasicMaterial color="#3b82f6" opacity={0.3} transparent />
        </mesh>
      )}
    </mesh>
  )
}
