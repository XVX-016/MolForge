
import { useEffect, useState, useMemo } from 'react';
import { useStudioStore } from '../../store/studioStore';
import { useHistoryStore } from '../../store/historyStore';
import Studio3DScene from './Studio3DScene';
import Card from '../ui/Card';
import { Plus, Trash2, Atom } from 'lucide-react';
import { addAtom, deleteAtom, deleteBond } from '../../lib/mutations';
import OptimizationPanel from './OptimizationPanel';
import SimulationTimeline from './SimulationTimeline';
import { generateVibrationTrajectory, type Trajectory, type SimulationFrame } from '../../lib/simulation';

export default function MolecularWorkspace() {
    const { mode, selection, setSelection } = useStudioStore();
    const { present: molecule, applyMutation } = useHistoryStore();

    const isEditable = mode === 'design';

    // --- Simulation State ---
    const [simulationTrajectory, setSimulationTrajectory] = useState<Trajectory | null>(null);
    const [simulatedFrame, setSimulatedFrame] = useState<SimulationFrame | null>(null);

    // Generate trajectory when entering simulate mode
    useEffect(() => {
        if (mode === 'simulate' && molecule) {
            // Generate a fresh trajectory
            const traj = generateVibrationTrajectory(molecule);
            setSimulationTrajectory(traj);
            setSimulatedFrame(traj[0]);
        } else {
            // Cleanup when leaving
            setSimulationTrajectory(null);
            setSimulatedFrame(null);
        }
    }, [mode, molecule]);


    // Which molecule to render?
    // In Simulate mode: use the frame's graph
    // Otherwise: use the store's molecule
    const moleculeToRender = (mode === 'simulate' && simulatedFrame)
        ? simulatedFrame.graph
        : molecule;


    // --- Actions ---

    const handleAddAtom = (element: string) => {
        if (!molecule) return;

        // Default position or offset from selection
        let position: [number, number, number] = [0, 0, 0];

        // Simple offset logic: if an atom is selected, place new atom nearby
        if (selection.type === 'atom' && selection.id) {
            const selectedAtom = molecule.atoms.find(a => a.id === selection.id);
            if (selectedAtom) {
                position = [
                    selectedAtom.position[0] + 1.5,
                    selectedAtom.position[1],
                    selectedAtom.position[2]
                ];
            }
        } else if (molecule.atoms.length > 0) {
            // If nothing selected but atoms exist, find a spot (naive)
            const lastAtom = molecule.atoms[molecule.atoms.length - 1];
            position = [
                lastAtom.position[0] + 1.5,
                lastAtom.position[1],
                lastAtom.position[2]
            ];
        }

        const newGraph = addAtom(molecule, element, position);
        applyMutation(newGraph, `Added ${element} Atom`, 'user');
    };

    const handleDelete = () => {
        if (!molecule || !selection.type || !selection.id) return;

        let newGraph = molecule;
        let description = '';

        if (selection.type === 'atom') {
            newGraph = deleteAtom(molecule, selection.id);
            description = 'Deleted Atom';
        } else if (selection.type === 'bond') {
            newGraph = deleteBond(molecule, selection.id);
            description = 'Deleted Bond';
        }

        applyMutation(newGraph, description, 'user');
        setSelection(null, null); // Clear selection
    };

    // --- Keyboard Shortcuts ---
    useEffect(() => {
        if (!isEditable) return;

        const handleKeyDown = (e: KeyboardEvent) => {
            if (e.key === 'Delete' || e.key === 'Backspace') {
                handleDelete();
            }
            // Deselect on Escape
            if (e.key === 'Escape') {
                setSelection(null, null);
            }
        };

        window.addEventListener('keydown', handleKeyDown);
        return () => window.removeEventListener('keydown', handleKeyDown);
    }, [isEditable, selection, molecule, applyMutation]);


    return (
        <Card className="h-full w-full overflow-hidden bg-gradient-to-br from-offwhite to-white border-lightGrey relative p-0 group">
            {/* 3D Scene Layer */}
            <div className="absolute inset-0 z-0 transition-opacity duration-300"
                style={{ pointerEvents: isEditable ? 'auto' : 'none' }}>
                <Studio3DScene
                    mode={mode}
                    molecule={moleculeToRender}
                    editable={isEditable}
                />
            </div>

            {/* Mode Status Overlay */}
            <div className="absolute top-4 left-4 z-10 pointer-events-none">
                <div className="inline-flex items-center gap-2 px-3 py-1.5 bg-white/90 backdrop-blur border border-lightGrey rounded-lg shadow-sm">
                    <div className={`w-2 h-2 rounded-full ${mode === 'design' ? 'bg-green-500 animate-pulse' :
                            mode === 'optimize' ? 'bg-purple-500' :
                                'bg-orange-500'
                        }`} />
                    <span className="text-xs font-medium text-darkGrey">
                        {mode === 'design' ? 'Live Editing Active' :
                            mode === 'optimize' ? 'Optimization Locked' :
                                'Simulation View'}
                    </span>
                </div>
            </div>

            {/* DESIGN TOOLBAR (Only in Design Mode) */}
            {mode === 'design' && (
                <div className="absolute top-4 right-4 z-20 flex flex-col gap-2">
                    <div className="bg-white rounded-xl shadow-lg border border-lightGrey p-2 flex flex-col gap-2">
                        <span className="text-[10px] font-bold text-center text-midGrey uppercase tracking-wider mb-1">Add</span>

                        <button
                            onClick={() => handleAddAtom('C')}
                            className="w-10 h-10 flex items-center justify-center rounded-lg bg-gray-50 hover:bg-black hover:text-white transition-colors border border-lightGrey font-bold text-sm"
                            title="Add Carbon"
                        >
                            C
                        </button>
                        <button
                            onClick={() => handleAddAtom('N')}
                            className="w-10 h-10 flex items-center justify-center rounded-lg bg-gray-50 hover:bg-black hover:text-white transition-colors border border-lightGrey font-bold text-sm"
                            title="Add Nitrogen"
                        >
                            N
                        </button>
                        <button
                            onClick={() => handleAddAtom('O')}
                            className="w-10 h-10 flex items-center justify-center rounded-lg bg-gray-50 hover:bg-black hover:text-white transition-colors border border-lightGrey font-bold text-sm"
                            title="Add Oxygen"
                        >
                            O
                        </button>
                        <div className="h-px bg-lightGrey w-full my-0.5" />
                        <button
                            onClick={handleDelete}
                            disabled={!selection.id}
                            className={`w-10 h-10 flex items-center justify-center rounded-lg transition-colors border border-lightGrey ${selection.id
                                    ? 'bg-red-50 text-red-600 hover:bg-red-600 hover:text-white cursor-pointer'
                                    : 'bg-transparent text-gray-300 cursor-not-allowed'
                                }`}
                            title="Delete Selected (Del)"
                        >
                            <Trash2 size={16} />
                        </button>
                    </div>

                    {selection.id && (
                        <div className="bg-black/80 text-white px-3 py-1.5 rounded-lg backdrop-blur text-xs font-mono text-center">
                            Selected: {selection.type}
                        </div>
                    )}
                </div>
            )}

            {/* OPTIMIZE MODE OVERLAY - Only if NOT in simulation */}
            {mode === 'optimize' && <OptimizationPanel />}

            {/* SIMULATE MODE TIMELINE */}
            {mode === 'simulate' && simulationTrajectory && (
                <SimulationTimeline
                    trajectory={simulationTrajectory}
                    onFrameChange={setSimulatedFrame}
                />
            )}
        </Card>
    );
}
