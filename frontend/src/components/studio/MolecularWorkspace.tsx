
import { useEffect, useState } from 'react';
import { useStudioStore } from '../../store/studioStore';
import { useHistoryStore } from '../../store/historyStore';
import { useStudioMode } from '../../lib/studio/hooks';
import { useStudioLayout } from '../../lib/studio/layout';
import Studio3DScene from './Studio3DScene';
import OptimizeCompareView from './OptimizeCompareView';
import { motion, AnimatePresence } from 'framer-motion';
import { deleteAtom, deleteBond } from '../../lib/mutations';
import OptimizationPanel from './OptimizationPanel';
import SimulationTimeline from './SimulationTimeline';
import { generateVibrationTrajectory, type Trajectory, type SimulationFrame } from '../../lib/simulation';

export default function MolecularWorkspace() {
    const { mode, selection, setSelection } = useStudioStore();
    const { present: molecule, applyMutation } = useHistoryStore();
    const { canEdit } = useStudioMode();
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
    // handleAddAtom removed as it's now orchestrated by AI

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
        <div className="h-full w-full overflow-hidden bg-[#F9FAFB] relative group">
            {/* 3D Scene Layer */}
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

            {/* Scientific Overlays (Logic-Bound) */}
            {layout.showOptimizationPanel && (
                <div className="absolute inset-x-0 bottom-12 z-20 flex justify-center pointer-events-none">
                    <div className="w-[500px] pointer-events-auto">
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

            {/* Selection HUD - Zero UI Aesthetic */}
            <AnimatePresence>
                {selection.id && (
                    <motion.div
                        initial={{ opacity: 0, scale: 0.95 }}
                        animate={{ opacity: 1, scale: 1 }}
                        exit={{ opacity: 0, scale: 0.95 }}
                        className="absolute bottom-12 left-1/2 -translate-x-1/2 px-5 py-2.5 bg-white/80 backdrop-blur-xl border border-gray-200 rounded-full flex items-center gap-4 pointer-events-none shadow-lg"
                    >
                        <div className="w-1.5 h-1.5 rounded-full bg-blue-500 shadow-[0_0_12px_rgba(59,130,246,0.5)]" />
                        <span className="text-[10px] font-black uppercase tracking-[0.2em] text-gray-700">
                            {selection.type === 'atom' ? `Atom ${selection.id}` : `Bond ${selection.id}`}
                        </span>
                    </motion.div>
                )}
            </AnimatePresence>
        </div>
    );
}
