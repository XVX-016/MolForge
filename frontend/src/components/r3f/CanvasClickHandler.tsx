import { useThree } from '@react-three/fiber'
import { useEffect } from 'react'
import { useMoleculeStore } from '../../store/moleculeStore'
import { addAtom } from '../../lib/engineAdapter'
import { screenToWorld } from '../../lib/raycasting'

interface PointerLikeMouseEvent {
  clientX: number
  clientY: number
  target: EventTarget | null
}

/**
 * Handles canvas clicks for add-atom tool
 */
export function CanvasClickHandler() {
  const { camera, gl } = useThree()
  const tool = useMoleculeStore((state) => state.tool)
  const atomToAdd = useMoleculeStore((state) => state.atomToAdd)

  useEffect(() => {
    if ((tool !== 'add-atom' && tool !== 'add_atom') || !atomToAdd) return

    const handleClick = (e: MouseEvent) => {
      // Only handle if clicking canvas (not atoms/bonds)
      if (e.target !== gl.domElement) return

      const canvas = gl.domElement

      // Convert to world coordinates
      const worldPos = screenToWorld(
        {
          target: canvas,
          clientX: e.clientX,
          clientY: e.clientY,
        } satisfies PointerLikeMouseEvent,
        camera,
        0
      )

      addAtom(atomToAdd, [worldPos.x, worldPos.y, worldPos.z])
    }

    gl.domElement.addEventListener('click', handleClick)
    return () => {
      gl.domElement.removeEventListener('click', handleClick)
    }
  }, [tool, atomToAdd, camera, gl])

  return null
}

