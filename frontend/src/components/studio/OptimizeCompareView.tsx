
import React from 'react';
import Studio3DScene from './Studio3DScene';
import type { MoleculeGraph } from '../../types/molecule';
import { useStudioStore } from '../../store/studioStore';
import { useHistoryStore } from '../../store/historyStore';
import { ArrowRight } from 'lucide-react';

export default function OptimizeCompareView() {
    const { mode } = useStudioStore();
    const { present: currentMolecule, past } = useHistoryStore();

    // The "Before" state is the first entry in history, or the current if none
    const beforeMolecule = past.length > 0 ? past[0] : currentMolecule;
    const afterMolecule = currentMolecule;

    if (!beforeMolecule || !afterMolecule) return null;

    return (
        <div className="flex w-full h-full gap-4 p-4 animate-in fade-in duration-500">
            {/* LEFT: ORIGINAL */}
            <div className="flex-1 relative rounded-3xl overflow-hidden border border-lightGrey bg-white shadow-sm">
                <div className="absolute top-4 left-6 z-10 px-3 py-1 bg-black/5 backdrop-blur-md rounded-lg">
                    <span className="text-[10px] font-black uppercase tracking-widest text-darkGrey">Original</span>
                </div>
                <Studio3DScene
                    mode={mode}
                    molecule={beforeMolecule}
                />
            </div>

            {/* DIVIDER: Transition Indicator */}
            <div className="flex flex-col items-center justify-center shrink-0 w-12 text-lightGrey/50">
                <ArrowRight size={24} className="animate-pulse" />
            </div>

            {/* RIGHT: CANDIDATE */}
            <div className="flex-1 relative rounded-3xl overflow-hidden border border-purple-200 bg-white shadow-xl shadow-purple-500/5 ring-4 ring-purple-500/5">
                <div className="absolute top-4 left-6 z-10 px-3 py-1 bg-purple-600 rounded-lg shadow-lg">
                    <span className="text-[10px] font-black uppercase tracking-widest text-white">Candidate</span>
                </div>
                <Studio3DScene
                    mode={mode}
                    molecule={afterMolecule}
                />
            </div>
        </div>
    );
}
