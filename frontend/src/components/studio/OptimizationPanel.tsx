
import React, { useMemo } from 'react';
import { useHistoryStore } from '../../store/historyStore';
import { generateSuggestions } from '../../lib/optimization';
import { Sparkles, ArrowRight, Zap } from 'lucide-react';

export default function OptimizationPanel() {
    const { present: molecule, applyMutation } = useHistoryStore();

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
        <div className="absolute top-4 right-4 bottom-4 w-80 flex flex-col gap-4 pointer-events-none">
            {/* Header */}
            <div className="bg-black/80 text-white p-4 rounded-xl backdrop-blur-md shadow-xl border border-white/10 pointer-events-auto">
                <div className="flex items-center gap-2 mb-2">
                    <Sparkles className="text-purple-400" size={18} />
                    <h3 className="font-bold text-sm uppercase tracking-wider">Optimization Engine</h3>
                </div>
                <p className="text-xs text-gray-300">
                    AI-driven suggestions to improve molecular properties.
                </p>
            </div>

            {/* Suggestion Cards */}
            <div className="flex-1 overflow-y-auto space-y-3 pr-2 scrollbar-hide pointer-events-auto">
                {suggestions.map((s) => (
                    <div key={s.id} className="bg-white/90 backdrop-blur border border-white shadow-lg rounded-xl p-4 hover:scale-[1.02] transition-transform duration-200 group">
                        <div className="flex justify-between items-start mb-2">
                            <h4 className="font-bold text-sm text-black">{s.title}</h4>
                            <span className="bg-purple-100 text-purple-700 text-[10px] font-bold px-1.5 py-0.5 rounded uppercase">
                                {s.type}
                            </span>
                        </div>

                        <p className="text-xs text-gray-600 mb-3 leading-relaxed">
                            {s.description}
                        </p>

                        {/* Impact Metrics */}
                        <div className="flex gap-3 mb-3 bg-gray-50 p-2 rounded-lg border border-gray-100">
                            <div className="flex flex-col">
                                <span className="text-[9px] text-gray-400 uppercase font-bold">LogP</span>
                                <span className={`text-xs font-mono font-bold ${s.impact.logPDelta > 0 ? 'text-orange-600' : 'text-blue-600'}`}>
                                    {s.impact.logPDelta > 0 ? '+' : ''}{s.impact.logPDelta}
                                </span>
                            </div>
                            <div className="flex flex-col">
                                <span className="text-[9px] text-gray-400 uppercase font-bold">Weight</span>
                                <span className="text-xs font-mono font-bold text-gray-700">
                                    +{s.impact.weightDelta.toFixed(1)}
                                </span>
                            </div>
                        </div>

                        <button
                            onClick={() => handleApply(s.id)}
                            className="w-full py-2 bg-black text-white text-xs font-bold rounded-lg flex items-center justify-center gap-2 group-hover:bg-purple-600 transition-colors"
                        >
                            <Zap size={12} /> Apply Optimization
                        </button>
                    </div>
                ))}
            </div>
        </div>
    );
}
