import React, { useRef, useEffect, useState } from 'react'
import * as THREE from 'three'
import gsap from 'gsap'
import type { Bond, Atom } from '../../../types/molecule'
import { useLabStore } from '../../../store/labStore'

interface BondMeshEnhancedProps {
  bond: Bond
  atom1: Atom
  atom2: Atom
  index: number
}

/**
 * Enhanced bond mesh with double/triple bond support and hover glow
 */
export default function BondMeshEnhanced({ bond, atom1, atom2, index }: BondMeshEnhancedProps) {
  const meshRef = useRef<THREE.Mesh>(null)
  const [hovered, setHovered] = useState(false)
  const selectedBondId = useLabStore(s => s.selectedBondId)
  const isSelected = bond.id === selectedBondId

  // Validate data
  if (!atom1 || !atom2 || !atom1.position || !atom2.position) {
    return null
  }

  useEffect(() => {
    if (meshRef.current) {
      // Grow animation
      gsap.fromTo(
        meshRef.current.scale,
        { y: 0 },
        {
          y: 1,
          duration: 0.25,
          delay: index * 0.02,
          ease: 'power2.out',
        }
      )
    }
  }, [index])

  const start = new THREE.Vector3(
    Number(atom1.position[0]) || 0,
    Number(atom1.position[1]) || 0,
    Number(atom1.position[2]) || 0
  )
  const end = new THREE.Vector3(
    Number(atom2.position[0]) || 0,
    Number(atom2.position[1]) || 0,
    Number(atom2.position[2]) || 0
  )
  const dir = new THREE.Vector3().subVectors(end, start)
  const len = dir.length()
  
  if (len === 0) return null
  
  const mid = start.clone().add(end).multiplyScalar(0.5)
  const quaternion = new THREE.Quaternion().setFromUnitVectors(
    new THREE.Vector3(0, 1, 0),
    dir.clone().normalize()
  )

  const baseRadius = 0.18
  const spacing = 0.08
  const bondOrder = bond.order || 1

  // Hover glow effect
  const glowIntensity = hovered || isSelected ? 1.3 : 1.0
  const color = hovered || isSelected ? 0x4676ff : 0x9aa0a6

  return (
    <group position={mid} quaternion={quaternion}>
      {/* Single bond or first cylinder of double/triple */}
      <mesh
        ref={bondOrder === 1 ? meshRef : undefined}
        onPointerOver={(e) => {
          e.stopPropagation()
          setHovered(true)
        }}
        onPointerOut={() => setHovered(false)}
        userData={{ bondId: bond.id, type: 'bond' }}
      >
        <cylinderGeometry args={[baseRadius, baseRadius, len, 12]} />
        <meshStandardMaterial
          color={color}
          metalness={0.0}
          roughness={0.6}
          emissive={hovered || isSelected ? 0x4676ff : 0x000000}
          emissiveIntensity={hovered || isSelected ? 0.2 : 0}
        />
      </mesh>

      {/* Second cylinder for double/triple bonds */}
      {bondOrder >= 2 && (
        <mesh position={[spacing, 0, 0]}>
          <cylinderGeometry args={[baseRadius * 0.9, baseRadius * 0.9, len, 12]} />
          <meshStandardMaterial
            color={color}
            metalness={0.0}
            roughness={0.6}
            emissive={hovered || isSelected ? 0x4676ff : 0x000000}
            emissiveIntensity={hovered || isSelected ? 0.2 : 0}
          />
        </mesh>
      )}

      {/* Third cylinder for triple bonds */}
      {bondOrder === 3 && (
        <mesh position={[-spacing, 0, 0]}>
          <cylinderGeometry args={[baseRadius * 0.9, baseRadius * 0.9, len, 12]} />
          <meshStandardMaterial
            color={color}
            metalness={0.0}
            roughness={0.6}
            emissive={hovered || isSelected ? 0x4676ff : 0x000000}
            emissiveIntensity={hovered || isSelected ? 0.2 : 0}
          />
        </mesh>
      )}
    </group>
  )
}

