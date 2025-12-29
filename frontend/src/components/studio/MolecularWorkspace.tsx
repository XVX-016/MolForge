
import { useEffect, useState } from 'react';
import { useStudioStore } from '../../store/studioStore';
import { useHistoryStore } from '../../store/historyStore';
import { useStudioMode } from '../../lib/studio/hooks';
import Studio3DScene from './Studio3DScene';
import Card from '../ui/Card';
import { Trash2, Atom, Info, Ruler, Box } from 'lucide-react';
import { addAtom, deleteAtom, deleteBond } from '../../lib/mutations';
import OptimizationPanel from './OptimizationPanel';
import SimulationTimeline from './SimulationTimeline';
import ToastOverlay from './ToastOverlay';
import { generateVibrationTrajectory, type Trajectory, type SimulationFrame } from '../../lib/simulation';

export default function MolecularWorkspace() {
    const { mode, selection, setSelection } = useStudioStore();
    const { present: molecule, applyMutation } = useHistoryStore();
    const { canEdit, modeColor } = useStudioMode();

    // --- Simulation State ---
    const [simulationTrajectory, setSimulationTrajectory] = useState<Trajectory | null>(null);
    const [simulatedFrame, setSimulatedFrame] = useState<SimulationFrame | null>(null);

    // Generate trajectory when entering simulate mode
    useEffect(() => {
        if (mode === 'simulate' && molecule) {
            const traj = generateVibrationTrajectory(molecule);
            setSimulationTrajectory(traj);
            setSimulatedFrame(traj[0]);
        } else {
            setSimulationTrajectory(null);
            setSimulatedFrame(null);
        }
    }, [mode, molecule]);

    const moleculeToRender = (mode === 'simulate' && simulatedFrame)
        ? simulatedFrame.graph
        : molecule;

    // --- Actions ---
    const handleAddAtom = (element: string) => {
        if (!molecule || !canEdit) return;

        let position: [number, number, number] = [0, 0, 0];
        if (selection.type === 'atom' && selection.id) {
            const selectedAtom = molecule.atoms.find(a => a.id === selection.id);
            if (selectedAtom) {
                position = [selectedAtom.position[0] + 1.5, selectedAtom.position[1], selectedAtom.position[2]];
            }
        } else if (molecule.atoms.length > 0) {
            const lastAtom = molecule.atoms[molecule.atoms.length - 1];
            position = [lastAtom.position[0] + 1.5, lastAtom.position[1], lastAtom.position[2]];
        }

        const newGraph = addAtom(molecule, element, position);
        applyMutation(newGraph, `Added ${element} Atom`, 'user');
    };

    const handleDelete = () => {
        if (!molecule || !selection.type || !selection.id || !canEdit) return;

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
        setSelection(null, null);
    };

    // --- Keyboard Shortcuts ---
    useEffect(() => {
        if (!canEdit) return;

        const handleKeyDown = (e: KeyboardEvent) => {
            if (e.key === 'Delete' || e.key === 'Backspace') {
                handleDelete();
            }
            if (e.key === 'Escape') {
                setSelection(null, null);
            }
        };

        window.addEventListener('keydown', handleKeyDown);
        return () => window.removeEventListener('keydown', handleKeyDown);
    }, [canEdit, selection, molecule, applyMutation]);

    return (
        <Card className="h-full w-full overflow-hidden bg-[#F8FAFC] border-lightGrey relative p-0 group">
            {/* 3D Scene Layer */}
            <div className="absolute inset-0 z-0">
                <Studio3DScene
                    mode={mode}
                    molecule={moleculeToRender}
                />
            </div>

            {/* Mode Status Overlay */}
            <div className="absolute top-6 left-6 z-10 pointer-events-none">
                <div className="inline-flex items-center gap-3 px-4 py-2.5 bg-white/95 backdrop-blur-xl border border-lightGrey/50 rounded-2xl shadow-xl">
                    <div className="w-2.5 h-2.5 rounded-full animate-pulse shadow-[0_0_10px_rgba(0,0,0,0.1)]" style={{ backgroundColor: modeColor }} />
                    <span className="text-[11px] font-black text-black uppercase tracking-[0.1em]">
                        {mode === 'design' ? 'Structure Authoring' :
                            mode === 'optimize' ? 'Optimization Phase' :
                                'Simulation Playback'}
                    </span>
                </div>
            </div>

            {/* Canvas Toolbar (Top Left, below badge) */}
            <div className="absolute top-20 left-6 z-10 flex gap-1.5 p-1.5 bg-white/80 backdrop-blur-md border border-lightGrey/30 rounded-2xl shadow-lg pointer-events-auto">
                <button className="p-2 rounded-xl bg-white shadow-sm border border-lightGrey/20 hover:bg-black hover:text-white transition-all text-midGrey" title="Ball & Stick View">
                    <Atom size={16} />
                </button>
                <div className="w-px h-6 bg-lightGrey/50 self-center mx-0.5" />
                <button className="p-2 rounded-xl hover:bg-gray-100 transition-all text-midGrey hover:text-black" title="Toggle Grid">
                    <Box size={16} />
                </button>
                <button className="p-2 rounded-xl hover:bg-gray-100 transition-all text-midGrey hover:text-black" title="Measurement Tool">
                    <Ruler size={16} />
                </button>
            </div>

            {/* DESIGN TOOLBAR (Only in Design Mode) */}
            {canEdit && (
                <div className="absolute top-6 right-6 z-20 flex flex-col gap-2 pointer-events-auto">
                    <div className="bg-white/95 backdrop-blur-xl rounded-2xl shadow-2xl border border-lightGrey/50 p-3 flex flex-col gap-2.5">
                        <span className="text-[9px] font-black text-center text-gray-400 uppercase tracking-[0.2em] mb-1">Editor</span>

                        <button
                            onClick={() => handleAddAtom('C')}
                            className="w-12 h-12 flex items-center justify-center rounded-xl bg-gray-50 hover:bg-black hover:text-white transition-all border border-lightGrey/30 font-black text-sm shadow-sm active:scale-90"
                            title="Add Carbon (C)"
                        >C</button>
                        <button
                            onClick={() => handleAddAtom('N')}
                            className="w-12 h-12 flex items-center justify-center rounded-xl bg-gray-50 hover:bg-black hover:text-white transition-all border border-lightGrey/30 font-black text-sm shadow-sm active:scale-90"
                            title="Add Nitrogen (N)"
                        >N</button>
                        <button
                            onClick={() => handleAddAtom('O')}
                            className="w-12 h-12 flex items-center justify-center rounded-xl bg-gray-50 hover:bg-black hover:text-white transition-all border border-lightGrey/30 font-black text-sm shadow-sm active:scale-90"
                            title="Add Oxygen (O)"
                        >O</button>

                        <div className="h-px bg-lightGrey/30 w-full my-1" />

                        <button
                            onClick={handleDelete}
                            disabled={!selection.id}
                            className={`w-12 h-12 flex items-center justify-center rounded-xl transition-all border border-lightGrey/30 shadow-sm active:scale-90 ${selection.id
                                ? 'bg-red-50 text-red-600 hover:bg-red-600 hover:text-white cursor-pointer'
                                : 'bg-gray-50/50 text-gray-200 cursor-not-allowed border-transparent shadow-none'
                                }`}
                            title="Delete Selected (Del)"
                        >
                            <Trash2 size={20} />
                        </button>
                    </div>
                </div>
            )}

            {/* Selection/Info Overlay (Bottom Right) */}
            {selection.id && (
                <div className="absolute bottom-6 right-6 z-10 animate-in fade-in slide-in-from-right-4 duration-300">
                    <div className="bg-black/90 text-white px-5 py-3 rounded-2xl backdrop-blur-xl border border-white/10 flex items-center gap-4 shadow-2xl">
                        <div className="p-1.5 bg-blue-500/20 rounded-lg">
                            <Info size={16} className="text-blue-400" />
                        </div>
                        <div className="flex flex-col">
                            <span className="text-[9px] uppercase tracking-[0.2em] font-black text-white/40 mb-0.5">Focus</span>
                            <span className="text-xs font-bold leading-none">{selection.type} <span className="text-blue-400">#{selection.id.slice(0, 6)}</span></span>
                        </div>
                    </div>
                </div>
            )}

            {/* Special Overlays per Mode */}
            {mode === 'optimize' && <OptimizationPanel />}
            {mode === 'simulate' && simulationTrajectory && (
                <SimulationTimeline
                    trajectory={simulationTrajectory}
                    onFrameChange={setSimulatedFrame}
                />
            )}

            {/* Feedback system */}
            <ToastOverlay />
        </Card>
    );
}
