
import { Suspense, useRef } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { Environment, Float, PerspectiveCamera, OrbitControls, Grid } from '@react-three/drei';
import * as THREE from 'three';
import type { MoleculeGraph } from '../../types/molecule';
import { moleculeToRenderable } from '../../lib/studioAdapter';
import type { StudioMode } from '../../types/studio';
import { useStudioStore } from '../../store/studioStore';
import { useStudioMode } from '../../lib/studio/hooks';
import type { SelectionType } from '../../store/studioStore';

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

// Internal component for the render logic
function FloatingMolecule({
    molecule,
    mode,
    onSelect,
    diffData
}: {
    molecule: MoleculeGraph,
    mode: StudioMode,
    onSelect?: (type: SelectionType, id: string | null) => void,
    diffData?: {
        atoms: { added: number[]; removed: number[]; modified: number[] };
        bonds: { added: number[][]; removed: number[][] };
    }
}) {
    const groupRef = useRef<THREE.Group>(null);
    const renderable = moleculeToRenderable(molecule);

    // Subscribe to store selection to show highlights
    const { selection } = useStudioStore();
    const { canEdit, canSelect, canOptimize } = useStudioMode();

    useFrame((state) => {
        if (!groupRef.current) return;
        const time = state.clock.elapsedTime;

        // Gentle float
        groupRef.current.position.y = Math.sin(time * 0.5) * 0.1;

        // Auto-rotate logic per mode
        if (mode === 'optimize') {
            groupRef.current.rotation.y += 0.001;
        }
    });

    return (
        <group ref={groupRef} scale={[0.5, 0.5, 0.5]}>
            {/* ATOMS */}
            {renderable.atoms.map((atom, i) => {
                const radius = ELEMENT_RADII[atom.element] || 1.0;

                // Diff Logic
                const isAdded = diffData?.atoms.added.includes(i);
                const isRemoved = diffData?.atoms.removed.includes(i);
                const isModified = diffData?.atoms.modified.includes(i);

                // Selection Logic
                const isSelected = canSelect && selection.type === 'atom' && selection.id === atom.id;

                // Optimization Highlight (Hover surrogate / Candidate atoms)
                // In a real version, this might come from the optimization engine's current focus
                const isCandidate = canOptimize && atom.element === 'C';

                // Visual Color Logic
                let color = "#ffffff";
                let emissive = "#000000";
                let emissiveIntensity = 0;

                if (isSelected) {
                    color = "#3b82f6"; // Blue selection
                    emissive = "#3b82f6";
                    emissiveIntensity = 0.5;
                } else if (isCandidate) {
                    color = "#a855f7"; // Purple hint
                    emissive = "#a855f7";
                    emissiveIntensity = 0.2;
                } else if (isAdded) {
                    color = "#22c55e"; // Green
                    emissive = "#22c55e";
                    emissiveIntensity = 0.4;
                } else if (isRemoved) {
                    color = "#ef4444"; // Red
                    emissive = "#ef4444";
                    emissiveIntensity = 0.4;
                } else if (isModified) {
                    color = "#f97316"; // Orange
                    emissive = "#f97316";
                    emissiveIntensity = 0.4;
                } else {
                    // Element Colors
                    if (atom.element === 'C') color = "#4B5563"; // Dark gray
                    if (atom.element === 'O') color = "#EF4444";
                    if (atom.element === 'N') color = "#2563EB";
                    if (atom.element === 'H') color = "#ffffff";
                }

                return (
                    <mesh
                        key={atom.id}
                        position={atom.position as any}
                        onClick={(e) => {
                            e.stopPropagation();
                            if ((canEdit || canOptimize) && onSelect) {
                                onSelect('atom', atom.id);
                            }
                        }}
                        onPointerOver={() => {
                            if (canEdit || canOptimize) document.body.style.cursor = 'pointer';
                        }}
                        onPointerOut={() => document.body.style.cursor = 'default'}
                    >
                        <sphereGeometry args={[radius, 32, 32]} />
                        <meshPhysicalMaterial
                            color={color}
                            metalness={0.4}
                            roughness={0.2}
                            envMapIntensity={1.5}
                            emissive={emissive}
                            emissiveIntensity={emissiveIntensity}
                        />
                        {/* Selection Halo */}
                        {isSelected && (
                            <mesh scale={[1.1, 1.1, 1.1]}>
                                <sphereGeometry args={[radius, 32, 32]} />
                                <meshBasicMaterial color="#60a5fa" transparent opacity={0.3} side={THREE.BackSide} />
                            </mesh>
                        )}
                    </mesh>
                );
            })}

            {/* BONDS */}
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

                const isSelected = canSelect && selection.type === 'bond' && selection.id === bond.id;

                // Diff Logic for bonds
                const a1Idx = parseInt(bond.fromAtomId.split('_')[1]);
                const a2Idx = parseInt(bond.toAtomId.split('_')[1]);

                const isBondAdded = diffData?.bonds.added.some(b =>
                    (b[0] === a1Idx && b[1] === a2Idx) || (b[0] === a2Idx && b[1] === a1Idx)
                );
                const isBondRemoved = diffData?.bonds.removed.some(b =>
                    (b[0] === a1Idx && b[1] === a2Idx) || (b[0] === a2Idx && b[1] === a1Idx)
                );

                let color = isSelected ? "#3b82f6" : "#E5E7EB";
                let emissive = isSelected ? "#3b82f6" : "#000000";
                const hasDiff = isBondAdded || isBondRemoved;

                if (isBondAdded) { color = "#22c55e"; emissive = "#22c55e"; }
                if (isBondRemoved) { color = "#ef4444"; emissive = "#ef4444"; }

                return (
                    <mesh
                        key={bond.id}
                        position={[mid.x, mid.y, mid.z]}
                        quaternion={[q.x, q.y, q.z, q.w]}
                        onClick={(e) => {
                            e.stopPropagation();
                            if (canEdit && onSelect) {
                                onSelect('bond', bond.id);
                            }
                        }}
                        onPointerOver={() => {
                            if (canEdit) document.body.style.cursor = 'pointer';
                        }}
                        onPointerOut={() => document.body.style.cursor = 'default'}
                    >
                        <cylinderGeometry args={[radius, radius, length, 32]} />
                        <meshPhysicalMaterial
                            color={color}
                            metalness={0.2}
                            roughness={0.1}
                            envMapIntensity={1}
                            emissive={emissive}
                            emissiveIntensity={(isSelected || hasDiff) ? 0.3 : 0}
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
    editable?: boolean;
    diffData?: {
        atoms: { added: number[]; removed: number[]; modified: number[] };
        bonds: { added: number[][]; removed: number[][] };
    };
    onCameraChange?: (target: THREE.Vector3, position: THREE.Vector3) => void;
    cameraTarget?: THREE.Vector3;
    cameraPosition?: THREE.Vector3;
}

export default function Studio3DScene({ molecule, mode, diffData, onCameraChange, cameraTarget, cameraPosition }: Studio3DSceneProps) {
    const { setSelection } = useStudioStore();
    const { canEdit, canOptimize } = useStudioMode();

    return (
        <div
            className={`w-full h-full relative ${(canEdit || canOptimize) ? 'cursor-move' : 'cursor-default'}`}
            onClick={() => {
                if (canEdit || canOptimize) setSelection(null, null);
            }}
        >
            <Canvas shadows dpr={[1, 2]}>
                <PerspectiveCamera
                    makeDefault
                    position={cameraPosition || [0, 0, 8]}
                    fov={45}
                />

                <Suspense fallback={null}>
                    <Environment preset="city" />
                    <ambientLight intensity={0.8} />
                    <pointLight position={[10, 10, 10]} intensity={1.5} castShadow />

                    {canEdit && (
                        <Grid infiniteGrid fadeDistance={30} sectionColor="#CBD5E1" cellColor="#E2E8F0" />
                    )}

                    <Float speed={1.5} rotationIntensity={0.2} floatIntensity={0.5}>
                        {molecule && (
                            <FloatingMolecule
                                molecule={molecule}
                                mode={mode}
                                onSelect={(type, id) => setSelection(type, id)}
                                diffData={diffData}
                            />
                        )}
                    </Float>

                    <OrbitControls
                        enableZoom={true}
                        enablePan={true}
                        autoRotate={mode === 'simulate'}
                        autoRotateSpeed={1.0}
                        enabled={true}
                        onChange={(e) => {
                            if (onCameraChange && e?.target) {
                                const controls = e.target as any;
                                onCameraChange(controls.target.clone(), controls.object.position.clone());
                            }
                        }}
                        target={cameraTarget}
                    />
                </Suspense>
            </Canvas>
        </div>
    );
}
