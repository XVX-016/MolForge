
import React, { useMemo } from 'react';
import { useHistoryStore } from '../../store/historyStore';
import { generateSuggestions } from '../../lib/optimization';
import { Sparkles, ArrowRight, Zap, Info } from 'lucide-react';
import { useStudioMode } from '../../lib/studio/hooks';

export default function OptimizationPanel() {
    const { present: molecule, applyMutation } = useHistoryStore();
    const { modeColor } = useStudioMode();

    const suggestions = useMemo(() => {
        if (!molecule) return [];
        return generateSuggestions(molecule);
    }, [molecule]);

    const handleApply = (strategyId: string) => {
        if (!molecule) return;
        const strategy = suggestions.find(s => s.id === strategyId);
        if (!strategy) return;

        const newGraph = strategy.apply(molecule);
        applyMutation(newGraph, `Applied Optimization: ${strategy.title}`, 'ai');
    };

    if (!molecule) return null;

    return (
        <div className="absolute top-6 right-6 bottom-6 w-80 flex flex-col gap-4 pointer-events-none animate-in fade-in slide-in-from-right-4 duration-500">
            {/* Header */}
            <div className="bg-black/90 text-white p-5 rounded-2xl backdrop-blur-xl shadow-2xl border border-white/10 pointer-events-auto">
                <div className="flex items-center gap-2 mb-2">
                    <Sparkles className="text-purple-400" style={{ color: modeColor }} size={16} />
                    <h3 className="font-black text-[10px] uppercase tracking-[0.2em]">Decision Engine</h3>
                </div>
                <p className="text-[11px] text-gray-400 leading-relaxed font-medium">
                    Property-steering suggestions based on active molecular geometry and valence constraints.
                </p>
            </div>

            {/* Suggestion Cards */}
            <div className="flex-1 overflow-y-auto space-y-3 pr-2 scrollbar-hide pointer-events-auto">
                {suggestions.map((s) => (
                    <div key={s.id} className="bg-white/95 backdrop-blur-md border border-lightGrey/50 shadow-xl rounded-2xl p-5 hover:scale-[1.02] transition-all duration-300 group ring-0 hover:ring-2 ring-purple-500/20">
                        <div className="flex justify-between items-start mb-2">
                            <h4 className="font-bold text-xs text-black uppercase tracking-tight">{s.title}</h4>
                            <span className="bg-purple-50 text-purple-700 text-[9px] font-black px-2 py-0.5 rounded-full uppercase tracking-widest border border-purple-100">
                                {s.type}
                            </span>
                        </div>

                        <p className="text-[11px] text-gray-500 mb-4 leading-relaxed font-medium">
                            {s.description}
                        </p>

                        {/* Rationale Section */}
                        <div className="mb-4 flex gap-2 items-start p-2.5 bg-gray-50 rounded-lg border border-gray-100">
                            <Info size={12} className="text-gray-400 mt-0.5" />
                            <p className="text-[10px] text-gray-400 italic font-medium">
                                Targeting valence saturation to stabilize steric strain.
                            </p>
                        </div>

                        {/* Impact Metrics */}
                        <div className="grid grid-cols-2 gap-3 mb-4">
                            <div className="flex flex-col p-2 bg-white border border-lightGrey/50 rounded-xl shadow-sm">
                                <span className="text-[9px] text-gray-400 uppercase font-black tracking-tighter">Δ LogP</span>
                                <span className={`text-xs font-mono font-bold ${s.impact.logPDelta > 0 ? 'text-orange-600' : 'text-blue-600'}`}>
                                    {s.impact.logPDelta > 0 ? '+' : ''}{s.impact.logPDelta}
                                </span>
                            </div>
                            <div className="flex flex-col p-2 bg-white border border-lightGrey/50 rounded-xl shadow-sm">
                                <span className="text-[9px] text-gray-400 uppercase font-black tracking-tighter">Δ Mass</span>
                                <span className="text-xs font-mono font-bold text-gray-900">
                                    +{s.impact.weightDelta.toFixed(1)}
                                </span>
                            </div>
                        </div>

                        <button
                            onClick={() => handleApply(s.id)}
                            className="w-full py-3 bg-black text-white text-[10px] font-black uppercase tracking-widest rounded-xl flex items-center justify-center gap-2 hover:bg-purple-600 transition-all shadow-lg active:scale-95"
                        >
                            <Zap size={12} /> Apply Mutation
                        </button>
                    </div>
                ))}
            </div>
        </div>
    );
}
