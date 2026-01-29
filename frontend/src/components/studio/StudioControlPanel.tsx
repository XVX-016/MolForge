
import { useState } from 'react';
import { Send, Sparkles, Zap, History, Info, ChevronRight } from 'lucide-react';
import { useStudioStore } from '../../store/studioStore';

export default function IntentPanel() {
    const [intent, setIntent] = useState('');
    const { status, dashboard, runAICommand, applyRule } = useStudioStore();

    const handleAISubmit = (e: React.FormEvent) => {
        e.preventDefault();
        if (!intent.trim() || status !== 'READY') return;
        runAICommand(intent);
        setIntent('');
    };

    return (
        <div className="flex flex-col h-full overflow-hidden bg-white">
            {/* Header: AI Orchestrator */}
            <div className="p-4 border-b border-[#E5E7EB] bg-slate-50/50">
                <div className="flex items-center justify-between mb-3">
                    <h3 className="text-[10px] font-black text-gray-500 uppercase tracking-widest flex items-center gap-2">
                        <Sparkles size={12} className="text-blue-600" /> AI Orchestrator
                    </h3>
                    <div className="flex items-center gap-1">
                        <div className={`w-1.5 h-1.5 rounded-full ${status === 'OPTIMIZING' ? 'bg-amber-500 animate-pulse' : 'bg-green-500'}`} />
                        <span className="text-[8px] font-bold text-gray-400">GEMINI 2.0</span>
                    </div>
                </div>

                <form onSubmit={handleAISubmit} className="relative">
                    <input
                        type="text"
                        value={intent}
                        onChange={(e) => setIntent(e.target.value)}
                        disabled={status !== 'READY'}
                        placeholder="Enter optimization intent..."
                        className="w-full bg-white border border-[#E5E7EB] rounded-lg pl-4 pr-10 py-3 text-[11px] outline-none focus:border-blue-500 focus:ring-4 focus:ring-blue-500/5 placeholder:text-gray-400 transition-all text-black font-medium"
                    />
                    <button
                        type="submit"
                        disabled={!intent.trim() || status !== 'READY'}
                        className="absolute right-2 top-1/2 -translate-y-1/2 p-2 rounded-md text-gray-400 hover:text-blue-600 disabled:opacity-30 transition-all hover:bg-blue-50"
                    >
                        <Send size={14} />
                    </button>
                </form>
                <p className="mt-2 text-[9px] text-gray-400 leading-tight">
                    AI will select deterministic rules. It cannot invent chemistry.
                </p>
            </div>

            {/* Rule List: Deterministic Tools */}
            <div className="flex-1 overflow-y-auto">
                <div className="p-4 border-b border-[#E5E7EB] bg-white sticky top-0 z-10">
                    <h4 className="text-[9px] font-black text-gray-400 uppercase tracking-widest flex items-center gap-2">
                        <Zap size={10} className="text-amber-500" /> Deterministic Rule Registry
                    </h4>
                </div>

                <div className="p-2 space-y-1">
                    {dashboard?.optimization_context.available_rules.map((rule) => (
                        <button
                            key={rule.id}
                            onClick={() => applyRule(rule.id)}
                            disabled={status !== 'READY'}
                            className="w-full text-left p-3 rounded-lg border border-transparent hover:border-[#E5E7EB] hover:bg-slate-50 transition-all group disabled:opacity-50 relative overflow-hidden"
                        >
                            <div className="flex justify-between items-start mb-1">
                                <span className="text-[11px] font-bold text-gray-900 group-hover:text-blue-600 transition-colors">{rule.title}</span>
                                <div className="flex items-center gap-1.5">
                                    {Object.entries(rule.impact).map(([key, val]) => (
                                        <div key={key} className="flex flex-col items-end">
                                            <span className={`text-[9px] font-mono leading-none ${val < 0 ? 'text-green-600' : 'text-blue-600'}`}>
                                                {val > 0 ? '+' : ''}{val}
                                            </span>
                                            <span className="text-[6px] uppercase font-black text-gray-300 tracking-tighter">{key.replace('Delta', '')}</span>
                                        </div>
                                    ))}
                                    <ChevronRight size={12} className="text-gray-300 group-hover:translate-x-0.5 transition-transform" />
                                </div>
                            </div>
                            <p className="text-[10px] text-gray-500 leading-tight mb-2 line-clamp-2">{rule.description}</p>
                            <div className="flex items-center gap-1.5 text-[8px] font-bold text-gray-400 uppercase">
                                <Info size={10} className="opacity-50" /> {rule.rationale}
                            </div>
                        </button>
                    ))}

                    {(!dashboard || dashboard.optimization_context.available_rules.length === 0) && (
                        <div className="flex flex-col items-center justify-center py-12 px-6 text-center">
                            <div className="w-8 h-8 rounded-full bg-slate-50 flex items-center justify-center mb-3">
                                <Zap size={14} className="text-slate-300" />
                            </div>
                            <p className="text-[10px] text-slate-400 font-medium italic">No applicable rules for this structure</p>
                        </div>
                    )}
                </div>
            </div>

            {/* Execution Timeline: Provenance */}
            <div className="p-4 border-t border-[#E5E7EB] bg-[#F9FAFB] shrink-0">
                <div className="flex items-center justify-between mb-4">
                    <h4 className="text-[9px] font-black text-gray-400 uppercase tracking-widest flex items-center gap-2">
                        <History size={10} /> Execution History
                    </h4>
                    <span className="text-[8px] font-bold text-blue-600 px-1.5 py-0.5 bg-blue-50 rounded uppercase">Audit Log Active</span>
                </div>

                <div className="space-y-3 relative before:absolute before:left-[5px] before:top-2 before:bottom-2 before:w-[1px] before:bg-slate-200">
                    <div className="flex items-start gap-3 relative">
                        <div className="w-2.5 h-2.5 rounded-full border-2 border-white bg-green-500 ring-4 ring-green-50 z-10 shrink-0 mt-0.5" />
                        <div className="flex-1 min-w-0">
                            <div className="flex justify-between items-baseline mb-0.5">
                                <span className="text-[10px] font-bold text-gray-900">Baseline Loaded</span>
                                <span className="text-[9px] font-mono text-gray-400 uppercase">14:22:01</span>
                            </div>
                            <p className="text-[9px] text-gray-500 truncate">Initial structural snapshot</p>
                        </div>
                    </div>

                    {status === 'OPTIMIZING' && (
                        <div className="flex items-start gap-3 relative animate-pulse">
                            <div className="w-2.5 h-2.5 rounded-full border-2 border-white bg-amber-500 ring-4 ring-amber-50 z-10 shrink-0 mt-0.5" />
                            <div className="flex-1 min-w-0">
                                <div className="flex justify-between items-baseline mb-0.5">
                                    <span className="text-[10px] font-bold text-amber-600">AI Orhcestration</span>
                                    <span className="text-[9px] font-mono text-amber-400 uppercase italic">Active...</span>
                                </div>
                                <p className="text-[9px] text-amber-500 truncate">Selecting rule from registry</p>
                            </div>
                        </div>
                    )}

                    {status === 'PROPOSED' && (
                        <div className="flex items-start gap-3 relative">
                            <div className="w-2.5 h-2.5 rounded-full border-2 border-white bg-blue-500 ring-4 ring-blue-50 z-10 shrink-0 mt-0.5" />
                            <div className="flex-1 min-w-0">
                                <div className="flex justify-between items-baseline mb-0.5">
                                    <span className="text-[10px] font-bold text-blue-600">Proposal Generated</span>
                                    <span className="text-[9px] font-mono text-blue-400 uppercase">Just Now</span>
                                </div>
                                <p className="text-[9px] text-blue-500 font-medium">Applied {dashboard?.optimization_context.available_rules[0]?.title}</p>
                            </div>
                        </div>
                    )}
                </div>
            </div>
        </div>
    );
}
