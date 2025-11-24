import React, { useState } from 'react'
import { Grid } from '@react-three/drei'
import { useLabStore } from '../../../store/labStore'

/**
 * Toggleable 3D grid helper
 */
export default function GridHelper() {
  const showGrid = useLabStore(s => s.showGrid ?? true)
  const setShowGrid = useLabStore(s => s.setShowGrid)

  if (!showGrid) return null

  return (
    <Grid
      args={[20, 20]}
      cellColor="#e0e0e0"
      sectionColor="#d0d0d0"
      cellThickness={0.5}
      sectionThickness={1}
      fadeDistance={25}
      fadeStrength={1}
      followCamera={false}
      infiniteGrid={true}
    />
  )
}

