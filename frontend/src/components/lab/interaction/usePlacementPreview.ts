import { useThree } from '@react-three/fiber'
import { useRef, useEffect, useState } from 'react'
import * as THREE from 'three'
import { useLabStore } from '../../../store/labStore'
import elementPalette from '../../../data/elementPalette'

/**
 * Hook for ghost preview during atom placement
 */
export function usePlacementPreview() {
  const { camera, raycaster, gl } = useThree()
  const [previewPos, setPreviewPos] = useState<THREE.Vector3 | null>(null)
  const currentTool = useLabStore(s => s.currentTool)
  const currentElement = useLabStore(s => s.currentElement)

  useEffect(() => {
    if (currentTool !== 'add_atom') {
      setPreviewPos(null)
      return
    }

    const handleMouseMove = (event: MouseEvent) => {
      const rect = gl.domElement.getBoundingClientRect()
      const pointer = new THREE.Vector2()
      pointer.x = ((event.clientX - rect.left) / rect.width) * 2 - 1
      pointer.y = -((event.clientY - rect.top) / rect.height) * 2 + 1

      // Raycast to y=0 plane
      raycaster.setFromCamera(pointer, camera)
      const plane = new THREE.Plane(new THREE.Vector3(0, 1, 0), 0)
      const point = new THREE.Vector3()
      raycaster.ray.intersectPlane(plane, point)

      if (point) {
        // Grid snapping (0.25 units)
        const snap = 0.25
        point.x = Math.round(point.x / snap) * snap
        point.y = Math.round(point.y / snap) * snap
        point.z = Math.round(point.z / snap) * snap
        setPreviewPos(point)
      }
    }

    gl.domElement.addEventListener('mousemove', handleMouseMove)
    return () => {
      gl.domElement.removeEventListener('mousemove', handleMouseMove)
    }
  }, [currentTool, camera, raycaster, gl])

  const palette = elementPalette[currentElement] || elementPalette.C
  const cpkColors: Record<string, number> = {
    H: 0xffffff, C: 0x909090, N: 0x3050f8, O: 0xff0d0d,
    F: 0x90e050, P: 0xff8000, S: 0xffff30, Cl: 0x1ff01f,
    Br: 0xa62929, I: 0x940094,
  }
  const color = new THREE.Color(cpkColors[currentElement] || 0xcccccc)
  const radius = palette?.radius || 1.0

  return { previewPos, color, radius }
}

