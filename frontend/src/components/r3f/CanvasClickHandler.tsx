import { useThree } from '@react-three/fiber'
import { useEffect } from 'react'
import { useMoleculeStore } from '../../store/moleculeStore'
import { addAtom } from '../../lib/engineAdapter'
import { screenToWorld } from '../../lib/raycasting'

/**
 * Handles canvas clicks for add-atom tool
 */
export function CanvasClickHandler() {
  const { camera, gl } = useThree()
  const tool = useMoleculeStore((state) => state.tool)
  const atomToAdd = useMoleculeStore((state) => state.atomToAdd)

  useEffect(() => {
    if (tool !== 'add-atom' || !atomToAdd) return

    const handleClick = (e: MouseEvent) => {
      // Only handle if clicking canvas (not atoms/bonds)
      if (e.target !== gl.domElement) return

      const canvas = gl.domElement
      const rect = canvas.getBoundingClientRect()

      // Convert to world coordinates
      const worldPos = screenToWorld(
        {
          ...e,
          target: canvas,
          clientX: e.clientX,
          clientY: e.clientY,
        } as any,
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

