
import { useState } from 'react';
import { Send, Sparkles, Zap, History, Info } from 'lucide-react';
import { useStudioStore } from '../../store/studioStore';

export default function StudioControlPanel() {
    const [intent, setIntent] = useState('');
    const { status, dashboard, runAICommand, applyRule } = useStudioStore();

    const handleAISubmit = (e: React.FormEvent) => {
        e.preventDefault();
        if (!intent.trim() || status !== 'READY') return;
        runAICommand(intent);
        setIntent('');
    };

    return (
        <div className="flex flex-col h-full overflow-hidden">
            <div className="p-4 border-b border-[#E5E7EB]">
                <h3 className="text-[10px] font-black text-gray-400 uppercase tracking-widest mb-3 flex items-center gap-2">
                    <Sparkles size={12} className="text-blue-600" /> Intent & Orchestration
                </h3>
                <form onSubmit={handleAISubmit} className="relative">
                    <input
                        type="text"
                        value={intent}
                        onChange={(e) => setIntent(e.target.value)}
                        disabled={status !== 'READY'}
                        placeholder="e.g. Reduce LogP but keep TPSA..."
                        className="w-full bg-white border border-[#E5E7EB] rounded-lg px-4 py-3 text-[11px] outline-none focus:border-blue-500 focus:ring-1 focus:ring-blue-500/10 placeholder:text-gray-400 transition-all text-black"
                    />
                    <button
                        type="submit"
                        disabled={!intent.trim() || status !== 'READY'}
                        className="absolute right-2 top-1/2 -translate-y-1/2 p-1.5 rounded-md text-gray-400 hover:text-blue-600 disabled:opacity-30 transition-all"
                    >
                        <Send size={14} />
                    </button>
                </form>
            </div>

            <div className="flex-1 overflow-y-auto p-4 space-y-4">
                <div>
                    <h4 className="text-[9px] font-black text-gray-400 uppercase tracking-wider mb-3 flex items-center gap-2">
                        <Zap size={10} className="text-amber-500" /> Available Transformations
                    </h4>
                    <div className="space-y-2">
                        {dashboard?.optimization_context.available_rules.map((rule) => (
                            <button
                                key={rule.id}
                                onClick={() => applyRule(rule.id)}
                                disabled={status !== 'READY'}
                                className="w-full text-left p-3 rounded-lg border border-[#E5E7EB] hover:border-blue-400 hover:bg-blue-50/30 transition-all group disabled:opacity-50"
                            >
                                <div className="flex justify-between items-start mb-1">
                                    <span className="text-[11px] font-bold text-gray-900">{rule.title}</span>
                                    <div className="flex items-center gap-1">
                                        {Object.entries(rule.impact).map(([key, val]) => (
                                            <span key={key} className={`text-[8px] font-mono ${val < 0 ? 'text-green-600' : 'text-blue-600'}`}>
                                                {key.slice(0, 4)}: {val > 0 ? '+' : ''}{val}
                                            </span>
                                        ))}
                                    </div>
                                </div>
                                <p className="text-[10px] text-gray-500 leading-tight mb-2">{rule.description}</p>
                                <div className="flex items-center gap-1.5 text-[8px] font-bold text-gray-400 uppercase group-hover:text-blue-500">
                                    <Info size={10} /> {rule.rationale.slice(0, 40)}...
                                </div>
                            </button>
                        ))}
                        {(!dashboard || dashboard.optimization_context.available_rules.length === 0) && (
                            <div className="text-center py-8 text-[11px] text-gray-400 italic">
                                No deterministic optimizations found for this structure.
                            </div>
                        )}
                    </div>
                </div>
            </div>

            <div className="p-4 border-t border-[#E5E7EB] bg-[#F9FAFB] shrink-0">
                <h4 className="text-[9px] font-black text-gray-400 uppercase tracking-wider mb-3 flex items-center gap-2">
                    <History size={10} /> Provenance Log
                </h4>
                <div className="space-y-1.5">
                    <div className="flex items-center justify-between text-[10px]">
                        <span className="text-gray-500">Loaded Baseline</span>
                        <span className="font-mono text-gray-400">14:22:01</span>
                    </div>
                    {status === 'PROPOSED' && (
                        <div className="flex items-center justify-between text-[10px] text-blue-600 font-bold">
                            <span>Rule Proposal Applied</span>
                            <span className="font-mono opacity-60 italic">Pending...</span>
                        </div>
                    )}
                </div>
            </div>
        </div>
    );
}
