import { useCallback, useEffect, useMemo, useState } from 'react'
import { Canvas, type ThreeEvent, useThree } from '@react-three/fiber'
import { Html, Line, OrbitControls } from '@react-three/drei'
import { useLabStore } from '../../store/labStore'
import AtomMesh from './AtomMesh'
import BondMesh from './BondMesh'

function CustomGrid() {
  const lines = useMemo(() => {
    const items = []
    const size = 100
    const spacing = 2
    const steps = size / spacing

    for (let i = -steps / 2; i <= steps / 2; i++) {
      const z = i * spacing
      items.push(
        <Line
          key={`h-${i}`}
          points={[[-size / 2, 0, z], [size / 2, 0, z]]}
          color="#e5e7eb"
          lineWidth={1}
          transparent
          opacity={0.4}
        />
      )
    }

    for (let i = -steps / 2; i <= steps / 2; i++) {
      const x = i * spacing
      items.push(
        <Line
          key={`v-${i}`}
          points={[[x, 0, -size / 2], [x, 0, size / 2]]}
          color="#e5e7eb"
          lineWidth={1}
          transparent
          opacity={0.4}
        />
      )
    }

    return items
  }, [])

  return <group position={[0, -0.01, 0]}>{lines}</group>
}

function Scene() {
  const { molecule, currentTool, currentElement, addAtom } = useLabStore()

  const onPlanePointerDown = useCallback(
    (event: ThreeEvent<PointerEvent>) => {
      if (currentTool !== 'add_atom') return

      event.stopPropagation()
      const point = event.point
      addAtom(currentElement, [point.x, point.y, point.z])
    },
    [addAtom, currentElement, currentTool]
  )

  return (
    <>
      <ambientLight intensity={0.9} />
      <directionalLight position={[5, 10, 7]} intensity={0.6} />
      <directionalLight position={[-5, 5, -5]} intensity={0.3} />

      <mesh rotation={[-Math.PI / 2, 0, 0]} position={[0, -0.01, 0]} receiveShadow onPointerDown={onPlanePointerDown}>
        <planeGeometry args={[100, 100]} />
        <meshBasicMaterial transparent opacity={0} />
      </mesh>

      <CustomGrid />

      <mesh position={[0, 0, 0]} rotation={[-Math.PI / 2, 0, 0]}>
        <circleGeometry args={[0.05, 32]} />
        <meshBasicMaterial color="#9ca3af" opacity={0.5} transparent />
      </mesh>

      {molecule.bonds.map((bond) => {
        const atomA = molecule.atoms.find((atom) => atom.id === bond.from)
        const atomB = molecule.atoms.find((atom) => atom.id === bond.to)

        if (!atomA || !atomB) {
          console.warn('[LabCanvas] Invalid bond or missing atom:', { bond, atomA, atomB })
          return null
        }

        const aPos = atomA.position
        const bPos = atomB.position

        return <BondMesh key={bond.id} aPos={[aPos.x, aPos.y, aPos.z]} bPos={[bPos.x, bPos.y, bPos.z]} order={bond.order} />
      })}

      {molecule.atoms.map((atom) => (
        <AtomMesh key={atom.id} atom={atom} />
      ))}
    </>
  )
}

function LabSceneWithControls({ onContextLost: onContextLostCallback }: { onContextLost?: () => void }) {
  const [contextLost, setContextLost] = useState(false)
  const { gl } = useThree()

  useEffect(() => {
    const canvas = gl.domElement as HTMLCanvasElement | undefined
    if (!canvas) return

    const onContextLost = (event: Event) => {
      event.preventDefault()
      console.warn('[LabCanvas] WebGL context lost')
      setContextLost(true)
      onContextLostCallback?.()
    }

    const onContextRestored = () => {
      console.info('[LabCanvas] WebGL context restored')
      setContextLost(false)
    }

    canvas.addEventListener('webglcontextlost', onContextLost)
    canvas.addEventListener('webglcontextrestored', onContextRestored)

    return () => {
      canvas.removeEventListener('webglcontextlost', onContextLost)
      canvas.removeEventListener('webglcontextrestored', onContextRestored)
    }
  }, [gl, onContextLostCallback])

  return (
    <>
      {contextLost && (
        <Html fullscreen>
          <div className="w-full h-full flex items-center justify-center bg-zinc-100 text-zinc-700">
            <div className="text-center space-y-2">
              <div className="font-semibold">Graphics context was lost</div>
              <div className="text-sm opacity-80">
                Try switching tabs more slowly or reloading the page. We&apos;ll attempt to recover automatically when possible.
              </div>
            </div>
          </div>
        </Html>
      )}

      {!contextLost && (
        <>
          <OrbitControls makeDefault enableDamping dampingFactor={0.1} maxPolarAngle={Math.PI / 2 - 0.1} minDistance={2} maxDistance={50} />
          <Scene />
        </>
      )}
    </>
  )
}

export default function LabCanvas() {
  const [canvasKey, setCanvasKey] = useState(0)
  const molecule = useLabStore((state) => state.molecule)

  const handleContextLost = useCallback(() => {
    console.warn('[LabCanvas] Forcing canvas remount due to context loss')
    setCanvasKey((prev) => prev + 1)
  }, [])

  return (
    <Canvas
      key={`canvas-${canvasKey}-${molecule.id}`}
      camera={{ position: [0, 6, 12], fov: 45 }}
      style={{ width: '100%', height: '100%' }}
      gl={{ alpha: true, antialias: true }}
      dpr={[1, 2]}
      frameloop="always"
    >
      <LabSceneWithControls onContextLost={handleContextLost} />
    </Canvas>
  )
}
