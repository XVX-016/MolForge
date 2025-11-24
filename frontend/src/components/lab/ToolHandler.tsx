import React from 'react'
import { useThree } from '@react-three/fiber'
import { useLabStore } from '../../store/labStore'
import { getTool } from '../../tools'
import { computeAutoBonds } from '../../utils/bondingEngine'
import * as THREE from 'three'
import { useCameraFocus } from './viewer/useCameraFocus'

/**
 * ToolHandler routes pointer events to the currently active tool
 */
export default function ToolHandler() {
  const { camera, gl, scene, raycaster } = useThree()
  const currentTool = useLabStore(s => s.currentTool)
  const store = useLabStore.getState()
  const autoBond = useLabStore(s => s.autoBond)
  const { focusOnAtom, resetCamera } = useCameraFocus()

  React.useEffect(() => {
    const handlePointerDown = (event: PointerEvent) => {
      // Only handle if clicking on canvas
      if (event.target !== gl.domElement && !gl.domElement.contains(event.target as Node)) return

      const rect = gl.domElement.getBoundingClientRect()
      const pointer = new THREE.Vector2()
      pointer.x = ((event.clientX - rect.left) / rect.width) * 2 - 1
      pointer.y = -((event.clientY - rect.top) / rect.height) * 2 + 1

      // Raycast to find clicked object
      raycaster.setFromCamera(pointer, camera)
      const intersects = raycaster.intersectObjects(scene.children, true)
      
      const r3fEvent = {
        clientX: event.clientX,
        clientY: event.clientY,
        target: gl.domElement,
        object: intersects.length > 0 ? intersects[0].object : null,
        camera: camera,
      }

      const tool = getTool(currentTool)
      if (tool && tool.onPointerDown) {
        tool.onPointerDown(r3fEvent, store)
        
        // Handle double-click to reset camera
        const now = Date.now()
        const lastClick = (window as any).__lastClick || 0
        if (now - lastClick < 300) {
          resetCamera()
        }
        ;(window as any).__lastClick = now

        // Focus on atom when selected (single click)
        if (currentTool === 'select' && r3fEvent.object?.userData?.atomId) {
          const atomId = r3fEvent.object.userData.atomId
          setTimeout(() => focusOnAtom(atomId), 100)
        }
        
        // Auto-bond after adding atom if enabled
        if (currentTool === 'add_atom' && autoBond) {
          setTimeout(() => {
            const mol = store.molecule
            const newBonds = computeAutoBonds(mol)
            newBonds.forEach(({ a, b }) => {
              store.addBond(a, b, 1)
            })
          }, 0)
        }
      }
    }

    gl.domElement.addEventListener('pointerdown', handlePointerDown)
    return () => {
      gl.domElement.removeEventListener('pointerdown', handlePointerDown)
    }
  }, [currentTool, autoBond, gl, camera, scene, raycaster, store])

  return null
}

