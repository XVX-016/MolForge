
import { Suspense, useRef, useEffect } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { Environment, Float, PerspectiveCamera, OrbitControls, Grid } from '@react-three/drei';
import * as THREE from 'three';
import type { MoleculeGraph } from '@biosynth/engine';
import { moleculeToRenderable } from '../../lib/engineAdapter';
import type { StudioMode } from '../../types/studio';

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
};

function FloatingMolecule({ molecule, mode }: { molecule: MoleculeGraph, mode: StudioMode }) {
    const groupRef = useRef<THREE.Group>(null);
    const renderable = moleculeToRenderable(molecule);

    useFrame((state) => {
        if (!groupRef.current) return;
        // Float effect - subtler in design mode
        const time = state.clock.elapsedTime;
        groupRef.current.position.y = Math.sin(time * 0.5) * 0.1;

        // Auto-rotate only in optimize/simulate modes as "presentation" view
        // In design mode, we want stability for editing (mocked)
        if (mode !== 'design') {
            groupRef.current.rotation.y += 0.002;
        }
    });

    return (
        <group ref={groupRef} scale={[0.5, 0.5, 0.5]}>
            {renderable.atoms.map((atom) => {
                const radius = ELEMENT_RADII[atom.element] || 1.0;
                // Highlight logic for Optimize Mode (Mock: Highlight Carbons)
                const isHighlight = mode === 'optimize' && atom.element === 'C';

                return (
                    <mesh key={atom.id} position={atom.position as any}>
                        <sphereGeometry args={[radius, 32, 32]} />
                        <meshPhysicalMaterial
                            color={isHighlight ? "#a855f7" : "#ffffff"} // Purple in optimize for Carbon
                            metalness={0.9}
                            roughness={0.1}
                            envMapIntensity={1.5}
                            emissive={isHighlight ? "#a855f7" : "#000000"}
                            emissiveIntensity={isHighlight ? 0.2 : 0}
                        />
                    </mesh>
                );
            })}

            {renderable.bonds.map((bond) => {
                const vFrom = new THREE.Vector3(...bond.from as [number, number, number]);
                const vTo = new THREE.Vector3(...bond.to as [number, number, number]);
                const diff = new THREE.Vector3().subVectors(vTo, vFrom);
                const length = diff.length();
                const mid = new THREE.Vector3().addVectors(vFrom, vTo).multiplyScalar(0.5);
                const q = new THREE.Quaternion().setFromUnitVectors(
                    new THREE.Vector3(0, 1, 0),
                    diff.clone().normalize()
                );
                const radius = bond.order === 1 ? 0.14 : bond.order === 2 ? 0.18 : 0.22;

                return (
                    <mesh
                        key={bond.id}
                        position={[mid.x, mid.y, mid.z]}
                        quaternion={[q.x, q.y, q.z, q.w]}
                    >
                        <cylinderGeometry args={[radius, radius, length, 32]} />
                        <meshPhysicalMaterial
                            color="#cccccc"
                            metalness={0.9}
                            roughness={0.1}
                            envMapIntensity={1.5}
                        />
                    </mesh>
                );
            })}
        </group>
    );
}

interface Studio3DSceneProps {
    mode: StudioMode;
    molecule: MoleculeGraph | null;
}

export default function Studio3DScene({ mode, molecule }: Studio3DSceneProps) {
    return (
        <div className="w-full h-full relative cursor-move">
            <Canvas shadows dpr={[1, 2]}>
                <PerspectiveCamera makeDefault position={[0, 0, 8]} fov={45} />

                <Suspense fallback={null}>
                    <Environment preset="city" />
                    <ambientLight intensity={0.5} />
                    <pointLight position={[10, 10, 10]} intensity={1} castShadow />

                    {/* Mode-specific Background/Grid */}
                    {mode === 'design' && (
                        <Grid infiniteGrid fadeDistance={30} sectionColor="#e5e7eb" cellColor="#f3f4f6" />
                    )}

                    <Float speed={1.5} rotationIntensity={0.2} floatIntensity={0.5}>
                        {molecule && (
                            <FloatingMolecule molecule={molecule} mode={mode} />
                        )}
                    </Float>

                    <OrbitControls
                        enableZoom={true}
                        enablePan={true}
                        autoRotate={false}
                    />
                </Suspense>
            </Canvas>
        </div>
    );
}
