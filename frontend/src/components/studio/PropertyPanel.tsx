
import React, { useMemo } from 'react';
import { useHistoryStore } from '../../store/historyStore';
import { calculateFormula, calculateMolecularWeight, estimateLogP, calculateTPSA } from '../../lib/chemistry';
import { Activity, Scale, Beaker, Zap } from 'lucide-react';
import Card from '../ui/Card';

export default function PropertyPanel() {
    // Subscribe to "present" (the current molecule)
    const { present: molecule } = useHistoryStore();

    // Deterministically recompute properties when molecule changes
    const stats = useMemo(() => {
        if (!molecule) return null;
        return {
            formula: calculateFormula(molecule),
            weight: calculateMolecularWeight(molecule),
            logP: estimateLogP(molecule),
            tpsa: calculateTPSA(molecule)
        };
    }, [molecule]);

    if (!stats) {
        return (
            <Card className="p-4 bg-white/50 border-lightGrey">
                <p className="text-xs text-midGrey text-center">No Molecule Loaded</p>
            </Card>
        );
    }

    return (
        <Card className="flex flex-col gap-4 bg-white border-lightGrey shadow-sm">
            <div className="flex items-center gap-2 border-b border-lightGrey pb-3">
                <Activity size={16} className="text-purple-600" />
                <span className="text-xs font-bold uppercase tracking-wider text-midGrey">
                    Properties (Live)
                </span>
            </div>

            <div className="grid grid-cols-2 gap-4">

                {/* Formula */}
                <div className="flex flex-col gap-1">
                    <span className="text-[10px] text-midGrey uppercase font-medium flex items-center gap-1">
                        <Beaker size={10} /> Formula
                    </span>
                    <span className="text-sm font-mono font-bold text-black tracking-tight">
                        {stats.formula}
                    </span>
                </div>

                {/* Weight */}
                <div className="flex flex-col gap-1">
                    <span className="text-[10px] text-midGrey uppercase font-medium flex items-center gap-1">
                        <Scale size={10} /> Weight
                    </span>
                    <span className="text-sm font-mono font-bold text-black tracking-tight">
                        {stats.weight} <span className="text-[10px] text-midGrey font-normal">g/mol</span>
                    </span>
                </div>

                {/* LogP */}
                <div className="flex flex-col gap-1">
                    <span className="text-[10px] text-midGrey uppercase font-medium flex items-center gap-1">
                        <Zap size={10} /> Est. LogP
                    </span>
                    <span className={`text-sm font-mono font-bold tracking-tight ${stats.logP > 5 ? 'text-orange-500' : 'text-black'}`}>
                        {stats.logP}
                    </span>
                </div>

                {/* TPSA */}
                <div className="flex flex-col gap-1">
                    <span className="text-[10px] text-midGrey uppercase font-medium flex items-center gap-1">
                        <Activity size={10} /> TPSA
                    </span>
                    <span className="text-sm font-mono font-bold text-black tracking-tight">
                        {stats.tpsa}
                    </span>
                </div>

            </div>
        </Card>
    );
}
