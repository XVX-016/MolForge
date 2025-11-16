/**
 * r3f-mock.tsx
 * Simple, local mock for three/react-three primitives that some components import directly.
 * (Some code imports './r3f-mock' or similar for test-time replacements.)
 */
import React from 'react';

export const Canvas = ({ children }: { children?: React.ReactNode }) => <div data-testid="canvas">{children}</div>;
export const mesh = (props: any) => <div data-testid="mesh" {...props} />;
export const sphereGeometry = (props: any) => <div data-testid="sphere-geometry" {...props} />;
export const meshStandardMaterial = (props: any) => <div data-testid="mesh-standard-material" {...props} />;
export const meshPhysicalMaterial = (props: any) => <div data-testid="mesh-physical-material" {...props} />;
export const ambientLight = (props: any) => <div data-testid="ambient-light" {...props} />;
export const directionalLight = (props: any) => <div data-testid="directional-light" {...props} />;
export const hemisphereLight = (props: any) => <div data-testid="hemisphere-light" {...props} />;
export const pointLight = (props: any) => <div data-testid="point-light" {...props} />;
export const cylinderGeometry = (props: any) => <div data-testid="cylinder-geometry" {...props} />;


