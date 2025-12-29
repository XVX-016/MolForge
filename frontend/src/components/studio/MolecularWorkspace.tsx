
import { useEffect, useState } from 'react';
import { useStudioStore } from '../../store/studioStore';
import { useHistoryStore } from '../../store/historyStore';
import { useStudioMode } from '../../lib/studio/hooks';
import { useStudioLayout } from '../../lib/studio/layout';
import Studio3DScene from './Studio3DScene';
import OptimizeCompareView from './OptimizeCompareView';
import Card from '../ui/Card';
import { Trash2, Info, Ruler, Box } from 'lucide-react';
import { addAtom, deleteAtom, deleteBond } from '../../lib/mutations';
import OptimizationPanel from './OptimizationPanel';
import SimulationTimeline from './SimulationTimeline';
import { generateVibrationTrajectory, type Trajectory, type SimulationFrame } from '../../lib/simulation';

export default function MolecularWorkspace() {
    const { mode, selection, setSelection } = useStudioStore();
    const { present: molecule, applyMutation } = useHistoryStore();
    const { canEdit, modeColor } = useStudioMode();
    const layout = useStudioLayout(mode);

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
            {/* 3D Scene Layer / Split View */}
            <div className="absolute inset-0 z-0">
                {layout.canvasLayout === 'split' ? (
                    <OptimizeCompareView />
                ) : (
                    <Studio3DScene
                        mode={mode}
                        molecule={moleculeToRender}
                    />
                )}
            </div>

            {/* Mode Status Badge (Differentiated) */}
            <div className="absolute top-8 left-8 z-10 pointer-events-none">
                <div className="inline-flex items-center gap-3 px-4 py-2 bg-white/90 backdrop-blur-xl border border-lightGrey/50 rounded-2xl shadow-xl">
                    <div className="w-2.5 h-2.5 rounded-full animate-pulse" style={{ backgroundColor: modeColor }} />
                    <span className="text-[10px] font-black text-black uppercase tracking-[0.2em]">
                        {mode === 'design' ? 'Authoring' :
                            mode === 'optimize' ? 'Optimization Console' :
                                'Simulation Playback'}
                    </span>
                </div>
            </div>

            {/* DESIGN TOOLBAR (Precision Authoring) - Only in Design Mode */}
            {layout.showEditorRail && (
                <div className="absolute top-8 right-8 z-20 flex flex-col gap-2 pointer-events-auto">
                    <div className="bg-white/95 backdrop-blur-xl rounded-2xl shadow-2xl border border-lightGrey/50 p-3 flex flex-col gap-2.5">
                        <span className="text-[9px] font-black text-center text-gray-400 uppercase tracking-[0.2em] mb-1">Elements</span>
                        <button onClick={() => handleAddAtom('C')} className="w-12 h-12 flex items-center justify-center rounded-xl bg-gray-50 hover:bg-black hover:text-white transition-all border border-lightGrey/30 font-black text-sm active:scale-90" title="Carbon">C</button>
                        <button onClick={() => handleAddAtom('N')} className="w-12 h-12 flex items-center justify-center rounded-xl bg-gray-50 hover:bg-black hover:text-white transition-all border border-lightGrey/30 font-black text-sm active:scale-90" title="Nitrogen">N</button>
                        <button onClick={() => handleAddAtom('O')} className="w-12 h-12 flex items-center justify-center rounded-xl bg-gray-50 hover:bg-black hover:text-white transition-all border border-lightGrey/30 font-black text-sm active:scale-90" title="Oxygen">O</button>
                        <div className="h-px bg-lightGrey/30 w-full my-1" />
                        <button onClick={handleDelete} disabled={!selection.id} className={`w-12 h-12 flex items-center justify-center rounded-xl transition-all border border-lightGrey/30 ${selection.id ? 'bg-red-50 text-red-600 hover:bg-red-600 hover:text-white' : 'bg-gray-50/50 text-gray-200 cursor-not-allowed border-transparent'}`} title="Delete Selection"><Trash2 size={20} /></button>
                    </div>

                    {/* View Controls below Editor Rail */}
                    <div className="bg-white/80 backdrop-blur-md rounded-2xl border border-lightGrey/30 p-2 flex flex-col gap-1.5 shadow-lg">
                        <button className="p-2.5 rounded-xl hover:bg-gray-100 text-midGrey transition-all" title="Measurement Tool"><Ruler size={16} /></button>
                        <button className="p-2.5 rounded-xl hover:bg-gray-100 text-midGrey transition-all" title="Toggle Grid"><Box size={16} /></button>
                    </div>
                </div>
            )}

            {/* Selection Info (Contextual) - Subtle Overlay */}
            {selection.id && mode !== 'simulate' && (
                <div className="absolute bottom-10 right-10 z-10 animate-in fade-in slide-in-from-right-4 duration-300">
                    <div className="bg-black/90 text-white px-5 py-3 rounded-2xl backdrop-blur-xl border border-white/10 flex items-center gap-4 shadow-2xl">
                        <div className="p-1.5 bg-blue-500/20 rounded-lg">
                            <Info size={16} className="text-blue-400" />
                        </div>
                        <div className="flex flex-col">
                            <span className="text-[9px] uppercase tracking-[0.2em] font-black text-white/40 mb-0.5">Selection</span>
                            <span className="text-xs font-bold leading-none">{selection.type} <span className="text-blue-400">#{selection.id.slice(0, 6)}</span></span>
                        </div>
                    </div>
                </div>
            )}

            {/* Mode-Specific Primary UI Overlays */}
            {layout.showOptimizationPanel && (
                <div className="absolute inset-x-0 bottom-0 z-20 p-8 flex justify-center pointer-events-none">
                    <div className="w-[600px] pointer-events-auto">
                        <OptimizationPanel />
                    </div>
                </div>
            )}

            {layout.showSimulationTimeline && simulationTrajectory && (
                <div className="absolute inset-x-0 bottom-0 z-20 p-8 bg-gradient-to-t from-black/20 to-transparent pointer-events-none">
                    <div className="pointer-events-auto w-full">
                        <SimulationTimeline
                            trajectory={simulationTrajectory}
                            onFrameChange={setSimulatedFrame}
                        />
                    </div>
                </div>
            )}
        </Card>
    );
}
