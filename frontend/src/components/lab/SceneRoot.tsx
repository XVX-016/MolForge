import React, { useEffect } from 'react'
import { Canvas, useThree } from '@react-three/fiber'
import { OrbitControls } from '@react-three/drei'
import gsap from 'gsap'
import MoleculeScene from './MoleculeScene'
import { useKeyboardShortcuts } from './hooks/useKeyboardShortcuts'
import GridHelper from './viewer/GridHelper'
import AtomGhostMesh from './viewer/AtomGhostMesh'

/**
 * Animates camera on mount
 * Must be used inside Canvas component
 */
function CameraAnimator() {
  try {
    const { camera } = useThree()

    useEffect(() => {
      // Store initial position
      const startPos = { x: 10, y: 10, z: 10 }
      camera.position.set(startPos.x, startPos.y, startPos.z)

      // Animate camera in
      gsap.to(camera.position, {
        x: 6,
        y: 6,
        z: 6,
        duration: 1.2,
        ease: 'power3.out',
      })
    }, [camera])

    return null
  } catch (error) {
    console.error('CameraAnimator must be used inside Canvas:', error)
    return null
  }
}

function KeyboardShortcuts() {
  useKeyboardShortcuts()
  return null
}

export default function SceneRoot(){
  return (
    <div style={{width:'100%', height:'100%', background: '#ffffff'}}>
      <KeyboardShortcuts />
      <Canvas 
        camera={{ position: [10, 10, 10], fov: 45 }}
        shadows
        style={{ background: '#ffffff' }}
        gl={{ 
          antialias: true,
          alpha: false,
          preserveDrawingBuffer: false,
          powerPreference: 'high-performance',
          stencil: false,
          depth: true
        }}
        dpr={[1, 2]}
        onCreated={({ gl }) => {
          gl.setPixelRatio(Math.min(window.devicePixelRatio, 2))
          // Handle WebGL context loss
          gl.domElement.addEventListener('webglcontextlost', (e) => {
            e.preventDefault()
            console.warn('WebGL context lost')
          })
          gl.domElement.addEventListener('webglcontextrestored', () => {
            console.log('WebGL context restored')
          })
        }}
      >
        <color attach="background" args={[0xffffff]} />
        <ambientLight intensity={0.7} />
        <directionalLight position={[10, 10, 5]} intensity={0.8} castShadow />
        <GridHelper />
        <React.Suspense fallback={null}>
          <MoleculeScene />
        </React.Suspense>
        <AtomGhostMesh />
        <CameraAnimator />
        <OrbitControls 
          ref={(ref) => {
            if (ref) {
              // Store ref for camera focus hook
              ;(window as any).__orbitControls = ref
            }
          }}
          makeDefault 
          enablePan={true}
          enableDamping
          dampingFactor={0.05}
          minDistance={2}
          maxDistance={50}
          rotateSpeed={0.5}
        />
      </Canvas>
    </div>
  )
}

