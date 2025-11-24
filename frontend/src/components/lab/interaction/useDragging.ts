import { useRef, useEffect } from 'react'
import { useThree } from '@react-three/fiber'
import * as THREE from 'three'
import { useLabStore } from '../../../store/labStore'

/**
 * Hook for dragging atoms in 3D space
 * Supports grid snapping and real-time bond updates
 */
export function useDragging() {
  const { camera, gl, raycaster } = useThree()
  const draggingRef = useRef<string | null>(null)
  const store = useLabStore.getState()

  useEffect(() => {
    const handlePointerMove = (event: PointerEvent) => {
      if (!draggingRef.current) return

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

        store.moveAtom(draggingRef.current, [point.x, point.y, point.z])
      }
    }

    const handlePointerUp = () => {
      draggingRef.current = null
      gl.domElement.style.cursor = 'default'
    }

    gl.domElement.addEventListener('pointermove', handlePointerMove)
    gl.domElement.addEventListener('pointerup', handlePointerUp)

    return () => {
      gl.domElement.removeEventListener('pointermove', handlePointerMove)
      gl.domElement.removeEventListener('pointerup', handlePointerUp)
    }
  }, [camera, gl, raycaster])

  const startDrag = (atomId: string) => {
    draggingRef.current = atomId
    gl.domElement.style.cursor = 'grabbing'
    store.snapshot() // Snapshot before drag starts
  }

  return { startDrag }
}

