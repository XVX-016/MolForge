import React from 'react'
import { useThree } from '@react-three/fiber'
import { useRef, useCallback } from 'react'
import * as THREE from 'three'
import gsap from 'gsap'
import { useLabStore } from '../../../store/labStore'

/**
 * Hook for camera focus and smooth transitions
 * Works with OrbitControls from drei
 */
export function useCameraFocus() {
  const { camera, scene } = useThree()
  const isAnimatingRef = useRef(false)
  const controlsRef = useRef<any>(null)

  // Get OrbitControls from window ref (set by SceneRoot)
  React.useEffect(() => {
    const checkControls = () => {
      if ((window as any).__orbitControls) {
        controlsRef.current = (window as any).__orbitControls
      }
    }
    checkControls()
    const interval = setInterval(checkControls, 100)
    return () => clearInterval(interval)
  }, [])

  const focusOnAtom = useCallback((atomId: string) => {
    if (isAnimatingRef.current) return

    const atom = useLabStore.getState().molecule.atoms.find(a => a.id === atomId)
    if (!atom) return

    const target = new THREE.Vector3(...atom.position)
    const distance = 8
    const offset = new THREE.Vector3(1, 1, 1).normalize().multiplyScalar(distance)
    const newPosition = target.clone().add(offset)

    isAnimatingRef.current = true

    gsap.to(camera.position, {
      x: newPosition.x,
      y: newPosition.y,
      z: newPosition.z,
      duration: 0.8,
      ease: 'power2.inOut',
      onComplete: () => {
        isAnimatingRef.current = false
      },
    })

    // Look at the atom
    if (controlsRef.current && controlsRef.current.target) {
      gsap.to(controlsRef.current.target, {
        x: target.x,
        y: target.y,
        z: target.z,
        duration: 0.8,
        ease: 'power2.inOut',
      })
    }
  }, [camera])

  const focusOnMolecule = useCallback(() => {
    if (isAnimatingRef.current) return

    const mol = useLabStore.getState().molecule
    if (mol.atoms.length === 0) return

    // Calculate bounding box
    const box = new THREE.Box3()
    mol.atoms.forEach(atom => {
      box.expandByPoint(new THREE.Vector3(...atom.position))
    })

    const center = box.getCenter(new THREE.Vector3())
    const size = box.getSize(new THREE.Vector3())
    const maxDim = Math.max(size.x, size.y, size.z)
    const distance = maxDim * 2.5

    const offset = new THREE.Vector3(1, 1, 1).normalize().multiplyScalar(distance)
    const newPosition = center.clone().add(offset)

    isAnimatingRef.current = true

    gsap.to(camera.position, {
      x: newPosition.x,
      y: newPosition.y,
      z: newPosition.z,
      duration: 1.0,
      ease: 'power2.inOut',
      onComplete: () => {
        isAnimatingRef.current = false
      },
    })

    if (controlsRef.current && controlsRef.current.target) {
      gsap.to(controlsRef.current.target, {
        x: center.x,
        y: center.y,
        z: center.z,
        duration: 1.0,
        ease: 'power2.inOut',
      })
    }
  }, [camera])

  const resetCamera = useCallback(() => {
    if (isAnimatingRef.current) return

    isAnimatingRef.current = true

    gsap.to(camera.position, {
      x: 10,
      y: 10,
      z: 10,
      duration: 1.0,
      ease: 'power2.inOut',
      onComplete: () => {
        isAnimatingRef.current = false
      },
    })

    if (controlsRef.current && controlsRef.current.target) {
      gsap.to(controlsRef.current.target, {
        x: 0,
        y: 0,
        z: 0,
        duration: 1.0,
        ease: 'power2.inOut',
      })
    }
  }, [camera])

  return { focusOnAtom, focusOnMolecule, resetCamera }
}
