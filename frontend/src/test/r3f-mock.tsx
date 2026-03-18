import React from 'react'

type PrimitiveProps = Record<string, unknown>

export const Canvas = ({ children }: { children?: React.ReactNode }) => <div data-testid="canvas">{children}</div>
export const mesh = (props: PrimitiveProps) => <div data-testid="mesh" {...props} />
export const sphereGeometry = (props: PrimitiveProps) => <div data-testid="sphere-geometry" {...props} />
export const meshStandardMaterial = (props: PrimitiveProps) => <div data-testid="mesh-standard-material" {...props} />
export const meshPhysicalMaterial = (props: PrimitiveProps) => <div data-testid="mesh-physical-material" {...props} />
export const ambientLight = (props: PrimitiveProps) => <div data-testid="ambient-light" {...props} />
export const directionalLight = (props: PrimitiveProps) => <div data-testid="directional-light" {...props} />
export const hemisphereLight = (props: PrimitiveProps) => <div data-testid="hemisphere-light" {...props} />
export const pointLight = (props: PrimitiveProps) => <div data-testid="point-light" {...props} />
export const cylinderGeometry = (props: PrimitiveProps) => <div data-testid="cylinder-geometry" {...props} />
