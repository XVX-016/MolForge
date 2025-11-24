import React from 'react'
import * as THREE from 'three'
import { usePlacementPreview } from '../interaction/usePlacementPreview'

/**
 * Ghost preview mesh for atom placement
 */
export default function AtomGhostMesh() {
  const { previewPos, color, radius } = usePlacementPreview()

  if (!previewPos) return null

  return (
    <mesh position={previewPos}>
      <sphereGeometry args={[radius, 32, 20]} />
      <meshStandardMaterial
        color={color}
        transparent
        opacity={0.25}
        metalness={0.0}
        roughness={0.6}
      />
    </mesh>
  )
}

