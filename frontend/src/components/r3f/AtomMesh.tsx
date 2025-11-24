import React, { useRef, useState } from 'react'
import * as THREE from 'three'
import { useFrame } from '@react-three/fiber'
import { Outlines, Text } from '@react-three/drei'
import { selectionManager } from './SelectionManager'
import { useMoleculeStore } from '../../store/moleculeStore'
import { kernelSelectionManager } from '../../kernel/selectionManager'

interface AtomMeshProps {
  id: string
  position: [number, number, number]
  element: 'C' | 'H' | 'O' | 'N' | 'F' | 'S' | 'P' | 'Cl' | 'Br' | 'I'
}

const ELEMENT_RADII: Record<string, number> = {
  H: 0.55,
  C: 1.0,
  O: 0.9,
  N: 0.95,
  F: 0.9,
  S: 1.2,
  P: 1.15,
  Cl: 1.1,
  Br: 1.25,
  I: 1.4,
}

// Professional metallic base colors (grey/ivory/chrome theme)
const ELEMENT_BASE_COLORS: Record<string, number> = {
  H: 0xF6F7F8,  // Soft white/ivory
  C: 0x9DA3AE,  // Charcoal grey
  O: 0x6B8FA3,  // Soft desaturated blue
  N: 0x8B95A1,  // Colder grey
  F: 0xA8B5B8,  // Soft beige-grey
  S: 0xC4B5A0,  // Soft beige tint
  P: 0xA8A39D,  // Warm grey
  Cl: 0x9FA8B3, // Light steel grey
  Br: 0x8B8B8B, // Medium grey
  I: 0x7A7A7A,  // Darker grey
}

// Subtle accent colors for element differentiation (thin ring/halo)
const ELEMENT_ACCENT_COLORS: Record<string, number> = {
  H: 0xFFFFFF,  // White edge
  C: 0x2B2E33,  // Dark charcoal ring
  O: 0x4A6FA5,  // Soft blue accent
  N: 0x5A6B7A,  // Cool grey accent
  F: 0x7A8B8F,  // Beige accent
  S: 0xB5A085,  // Beige accent
  P: 0x8B7A6B,  // Warm accent
  Cl: 0x6B7A8B, // Steel accent
  Br: 0x6B6B6B, // Grey accent
  I: 0x5A5A5A,  // Dark accent
}

export default function AtomMesh({ id, position, element }: AtomMeshProps) {
  const meshRef = useRef<THREE.Mesh>(null)
  const [isHovered, setIsHovered] = useState(false)
  const [isDragging, setIsDragging] = useState(false)
  const selectedAtomId = useMoleculeStore((state) => state.selectedAtomId)
  const tool = useMoleculeStore((state) => state.tool)
  
  const radius = ELEMENT_RADII[element] || 1.0
  const baseColor = ELEMENT_BASE_COLORS[element] || 0x9DA3AE
  const accentColor = ELEMENT_ACCENT_COLORS[element] || 0x2B2E33
  const isSelected = selectedAtomId === id
  
  // Visual feedback based on tool - using neonCyan for selection
  let scale = isSelected ? 1.15 : 1.0
  let outlineColor = 0x8BF3FF // neonCyan for select/hover
  if (tool === 'delete' && (isHovered || isSelected)) {
    outlineColor = 0xC6BDFE // violetEdge for delete
  } else if (tool === 'bond' && (isHovered || isSelected)) {
    outlineColor = 0x8BF3FF // neonCyan for bond
  }
  
  const opacity = isDragging ? 0.7 : 1.0

  // Handle pointer over (hover)
  const handlePointerOver = (e: React.PointerEvent) => {
    e.stopPropagation()
    setIsHovered(true)
    selectionManager.onHover(id)
  }

  // Handle pointer out (unhover)
  const handlePointerOut = (e: React.PointerEvent) => {
    e.stopPropagation()
    setIsHovered(false)
    if (!isDragging) {
      selectionManager.onHover(null)
    }
  }

  // Handle click (select or delete)
  const handleClick = (e: React.MouseEvent) => {
    e.stopPropagation()
    
    if (tool === 'delete') {
      // Import removeAtom dynamically to avoid circular dependency
      import('../../lib/engineAdapter').then(({ removeAtom: remove }) => {
        remove(id)
      })
    } else {
      // Update both UI and kernel selection managers
      selectionManager.onSelect(id)
      kernelSelectionManager.selectAtom(id)
      useMoleculeStore.getState().selectAtom(id)
    }
  }

  // Handle pointer down (start drag)
  const handlePointerDown = (e: React.PointerEvent) => {
    e.stopPropagation()
    setIsDragging(true)
    selectionManager.startDrag(id)
    // Prevent orbit controls
    ;(e.target as HTMLElement).setPointerCapture(e.pointerId)
  }

  // Handle pointer move (position updates handled by InteractionLayer)
  const handlePointerMove = (e: React.PointerEvent) => {
    if (!isDragging) return
    e.stopPropagation()
    // Position updates are handled by InteractionLayer component
    // which has access to camera and can properly convert coordinates
  }

  // Handle pointer up (end drag)
  const handlePointerUp = (e: React.PointerEvent) => {
    if (isDragging) {
      e.stopPropagation()
      setIsDragging(false)
      selectionManager.endDrag()
      ;(e.target as HTMLElement).releasePointerCapture(e.pointerId)
    }
  }

  // Update position from store
  useFrame(() => {
    if (meshRef.current) {
      meshRef.current.position.set(...position)
    }
  })

  const segments = (typeof window !== 'undefined' && window.devicePixelRatio && window.devicePixelRatio > 1.5) ? 48 : 64

  return (
    <mesh
      ref={meshRef}
      position={position}
      scale={scale}
      onClick={handleClick}
      onPointerOver={handlePointerOver}
      onPointerOut={handlePointerOut}
      onPointerDown={handlePointerDown}
      onPointerMove={handlePointerMove}
      onPointerUp={handlePointerUp}
      cursor={isDragging ? 'grabbing' : 'pointer'}
    >
      <sphereGeometry args={[radius, segments, segments]} />
      <meshPhysicalMaterial
        color={baseColor}
        roughness={0.2}
        metalness={0.4}
        clearcoat={1}
        clearcoatRoughness={0.1}
        thickness={0.5}
        transparent={isDragging}
        opacity={opacity}
        envMapIntensity={1.2}
      />
      {/* Subtle accent ring on hover/select */}
      {(isHovered || isSelected) && (
        <>
          <Outlines thickness={0.08} color={accentColor} />
          <Outlines thickness={0.12} color={outlineColor} />
        </>
      )}
      {/* Element symbol label - always visible, facing camera */}
      <Text
        position={[0, 0, radius * 1.05]}
        fontSize={radius * 0.5}
        color={isSelected || isHovered ? "#000000" : "#2B2E33"}
        anchorX="center"
        anchorY="middle"
        outlineWidth={0.01}
        outlineColor="#FFFFFF"
        renderOrder={1000}
      >
        {element}
      </Text>
    </mesh>
  )
}

