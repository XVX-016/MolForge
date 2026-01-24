import React from 'react';
import { useStudioStore } from '../../store/studioStore';
import { Sparkles, Zap, Info, AlertTriangle, CheckCircle } from 'lucide-react';
import { motion, AnimatePresence } from 'framer-motion';

export default function OptimizationPanel() {
    const { analysis, isComputingProperties, runCommand } = useStudioStore();

    if (!analysis) return (
        <div className="absolute top-6 right-6 bottom-6 w-80 flex flex-col justify-center items-center gap-4 pointer-events-none">
            <div className="bg-white/90 backdrop-blur-xl p-8 rounded-2xl border border-gray-100 shadow-sm text-center">
                <div className="w-10 h-10 bg-blue-50 rounded-full flex items-center justify-center mx-auto mb-4">
                    <Sparkles className="text-blue-500 w-5 h-5" />
                </div>
                <h3 className="text-sm font-bold text-gray-900 mb-1">Architecture Analysis</h3>
                <p className="text-xs text-gray-400">Add atoms to begin structural assessment.</p>
            </div>
        </div>
    );

    const { issues, suggestions } = analysis;

    const getSeverityColor = (severity: string) => {
        switch (severity) {
            case 'high': return 'text-red-600 bg-red-50 border-red-100';
            case 'medium': return 'text-amber-600 bg-amber-50 border-amber-100';
            default: return 'text-blue-600 bg-blue-50 border-blue-100';
        }
    };

    const getSeverityIcon = (severity: string) => {
        switch (severity) {
            case 'high': return <AlertTriangle size={14} />;
            case 'medium': return <Info size={14} />;
            default: return <CheckCircle size={14} />;
        }
    };

    return (
        <div className="absolute top-6 right-6 bottom-6 w-80 flex flex-col gap-4 pointer-events-none animate-in fade-in slide-in-from-right-4 duration-500 overflow-hidden font-sans">
            {/* Clinical Insights Header */}
            <div className="bg-gray-900 text-white p-5 rounded-2xl backdrop-blur-xl shadow-2xl border border-white/10 pointer-events-auto">
                <div className="flex items-center justify-between mb-2">
                    <div className="flex items-center gap-2">
                        <Sparkles className="text-blue-400" size={16} />
                        <h3 className="font-black text-[10px] uppercase tracking-[0.2em]">Clinical Insights</h3>
                    </div>
                    {isComputingProperties && (
                        <motion.div
                            animate={{ rotate: 360 }}
                            transition={{ repeat: Infinity, duration: 1, ease: "linear" }}
                            className="w-3 h-3 border border-blue-400 border-t-transparent rounded-full"
                        />
                    )}
                </div>
                <p className="text-[10px] text-gray-400 leading-relaxed font-medium">
                    Deterministic structural auditing via RDKit FilterCatalog and bioisosteric heuristic engine.
                </p>
            </div>

            <div className="flex-1 overflow-y-auto space-y-4 pr-2 scrollbar-hide pointer-events-auto pb-6">
                {/* Structural Alerts Section */}
                {issues.length > 0 && (
                    <div className="space-y-2">
                        <h4 className="text-[9px] font-black text-gray-400 uppercase tracking-[0.15em] ml-1">Structural Alerts</h4>
                        {issues.map((issue, idx) => (
                            <motion.div
                                initial={{ opacity: 0, scale: 0.95 }}
                                animate={{ opacity: 1, scale: 1 }}
                                key={idx}
                                className={`p-3 rounded-xl border flex gap-3 ${getSeverityColor(issue.severity)}`}
                            >
                                <div className="mt-0.5">{getSeverityIcon(issue.severity)}</div>
                                <div>
                                    <h5 className="font-bold text-[11px] leading-tight mb-0.5">{issue.title}</h5>
                                    <p className="text-[10px] opacity-80 leading-snug">{issue.description}</p>
                                </div>
                            </motion.div>
                        ))}
                    </div>
                )}

                {/* Optimization Suggestions Section */}
                {suggestions.length > 0 && (
                    <div className="space-y-3">
                        <h4 className="text-[9px] font-black text-gray-400 uppercase tracking-[0.15em] ml-1">Heuristic Steering</h4>
                        {suggestions.map((s) => (
                            <div key={s.id} className="bg-white/95 backdrop-blur-md border border-gray-200 shadow-sm rounded-xl p-4 hover:shadow-md transition-all group border-b-2 hover:border-blue-500">
                                <div className="flex justify-between items-start mb-2">
                                    <h4 className="font-bold text-xs text-gray-900 group-hover:text-blue-600 transition-colors uppercase tracking-tight">{s.title}</h4>
                                    <span className="text-[8px] font-black px-1.5 py-0.5 rounded bg-gray-100 text-gray-500 uppercase tracking-widest">
                                        Optimization
                                    </span>
                                </div>

                                <p className="text-[10px] text-gray-500 mb-3 leading-relaxed font-medium">
                                    {s.description}
                                </p>

                                {/* Property Impact Preview */}
                                <div className="grid grid-cols-2 gap-2 mb-3">
                                    <div className="flex flex-col p-1.5 bg-gray-50/50 rounded-lg border border-gray-100">
                                        <span className="text-[8px] text-gray-400 uppercase font-bold tracking-tighter">Δ LogP</span>
                                        <span className={`text-[10px] font-mono font-bold ${s.impact.logPDelta > 0 ? 'text-amber-600' : 'text-blue-600'}`}>
                                            {s.impact.logPDelta > 0 ? '+' : ''}{s.impact.logPDelta}
                                        </span>
                                    </div>
                                    <div className="flex flex-col p-1.5 bg-gray-50/50 rounded-lg border border-gray-100">
                                        <span className="text-[8px] text-gray-400 uppercase font-bold tracking-tighter">Δ TPSA</span>
                                        <span className={`text-[10px] font-mono font-bold ${s.impact.tpsaDelta > 0 ? 'text-emerald-600' : 'text-amber-600'}`}>
                                            {s.impact.tpsaDelta > 0 ? '+' : ''}{s.impact.tpsaDelta}
                                        </span>
                                    </div>
                                </div>

                                {/* AI Prompt Trigger */}
                                <button
                                    onClick={() => runCommand(`${s.title}: ${s.description}`)}
                                    className="w-full py-2 bg-gray-900 text-white text-[9px] font-black uppercase tracking-widest rounded-lg flex items-center justify-center gap-2 hover:bg-blue-600 transition-all active:scale-95"
                                >
                                    <Zap size={10} /> Propose Mutation
                                </button>
                            </div>
                        ))}
                    </div>
                )}

                {issues.length === 0 && suggestions.length === 0 && (
                    <div className="p-8 text-center bg-white/50 rounded-2xl border border-dashed border-gray-200">
                        <CheckCircle className="text-emerald-500 w-6 h-6 mx-auto mb-2 opacity-50" />
                        <p className="text-[10px] text-gray-400 font-medium">Structure satisfies all active heuristic filters.</p>
                    </div>
                )}
            </div>
        </div>
    );
}
