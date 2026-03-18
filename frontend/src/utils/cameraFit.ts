import { Box3, Vector3 } from 'three'
import type { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js'

interface AtomPositionLike {
  x: number
  y: number
  z?: number
}

interface PerspectiveCameraLike {
  fov: number
  aspect: number
  position: Vector3
  lookAt: (target: Vector3) => void
}

export function fitCameraToAtoms(
  camera: PerspectiveCameraLike,
  controls: OrbitControls | null | undefined,
  atomPositions: AtomPositionLike[],
  offset = 1.25
) {
  if (!atomPositions.length) return

  const box = new Box3()
  atomPositions.forEach((p) => box.expandByPoint(new Vector3(p.x, p.y, p.z ?? 0)))

  const size = box.getSize(new Vector3())
  const center = box.getCenter(new Vector3())
  const maxSize = Math.max(size.x, size.y, size.z)
  const fitHeightDistance = maxSize / (2 * Math.atan(Math.PI * camera.fov / 360))
  const fitWidthDistance = fitHeightDistance / camera.aspect
  const distance = offset * Math.max(fitHeightDistance, fitWidthDistance)
  const direction = new Vector3(0, 0, 1)

  camera.position.copy(direction.multiplyScalar(distance).add(center))
  camera.lookAt(center)

  if (controls) {
    controls.target.copy(center)
    controls.update()
  }
}
