import { Sparkles, FlaskConical, ShieldAlert, ChevronRight } from 'lucide-react';

interface GeminiDeskProps {
    insight?: {
        summary: string;
        key_observations: Array<{ type: string; impact: string; text: string }>;
        contradictions: Array<{ source1: string; source2: string; text: string }>;
        suggested_next_steps: Array<{ action: string; rationale: string }>;
    };
    loading?: boolean;
    isEmpty?: boolean;
}

export default function GeminiDesk({ insight, loading, isEmpty }: GeminiDeskProps) {
    if (isEmpty) {
        return (
            <div className="flex flex-col h-full bg-white">
                <div className="p-4 border-b border-[#E5E7EB] flex items-center justify-between">
                    <div>
                        <h3 className="text-xs font-black uppercase tracking-widest text-black flex items-center gap-2">
                            <Sparkles size={14} className="text-gray-400" />
                            Scientific Assistant
                        </h3>
                    </div>
                </div>
                <div className="p-8 flex flex-col items-center justify-center h-full text-center text-midGrey">
                    <p className="text-[11px] leading-relaxed max-w-[240px] text-gray-400 italic">
                        No analysis available yet.
                    </p>
                    <ul className="mt-4 text-[10px] text-left space-y-2 text-gray-500">
                        <li className="flex gap-2"><span>•</span> Interpret docking & MD results</li>
                        <li className="flex gap-2"><span>•</span> Flag stability or binding risks</li>
                        <li className="flex gap-2"><span>•</span> Suggest next scientific steps</li>
                    </ul>
                </div>
            </div>
        );
    }
    if (loading) {
        return (
            <div className="p-6 space-y-4 animate-pulse">
                <div className="h-4 w-1/3 bg-gray-100 rounded" />
                <div className="h-24 w-full bg-gray-100 rounded-2xl" />
                <div className="h-40 w-full bg-gray-100 rounded-2xl" />
            </div>
        );
    }

    if (!insight) {
        return (
            <div className="p-8 flex flex-col items-center justify-center h-full text-center text-midGrey">
                <div className="w-16 h-16 bg-gray-50 rounded-full flex items-center justify-center mb-4">
                    <Sparkles className="text-gray-200" size={32} />
                </div>
                <h4 className="text-xs font-black uppercase tracking-widest mb-2">Scientific Reasoning Desk</h4>
                <p className="text-[11px] leading-relaxed max-w-[240px]">
                    Select a completed analysis node to trigger AI interpretation and hypothesis generation.
                </p>
            </div>
        );
    }

    return (
        <div className="flex flex-col h-full bg-white">
            <div className="p-4 border-b border-[#E5E7EB] flex items-center justify-between">
                <div>
                    <h3 className="text-xs font-black uppercase tracking-widest text-black flex items-center gap-2">
                        <Sparkles size={14} className="text-blue-600" />
                        Scientific Interpretation
                    </h3>
                </div>
                <span className="text-[10px] font-bold text-blue-600 bg-blue-50 px-2 py-0.5 rounded-full">
                    Gemini 2.0
                </span>
            </div>

            <div className="flex-1 overflow-y-auto p-6 space-y-6">
                {/* Summary Section */}
                <section>
                    <p className="text-sm font-medium leading-relaxed text-black italic">
                        "{insight.summary}"
                    </p>
                </section>

                {/* Key Observations */}
                <section className="space-y-3">
                    <h4 className="text-[10px] font-black uppercase text-midGrey tracking-wider">Key Observations</h4>
                    {insight.key_observations.map((obs, idx) => (
                        <div key={idx} className="p-3 bg-[#F9FAFB] border border-[#E5E7EB] rounded-xl flex gap-3">
                            <FlaskConical size={16} className="text-blue-500 shrink-0" />
                            <div>
                                <span className={`text-[9px] font-bold uppercase ${obs.impact === 'high' ? 'text-red-600' : 'text-blue-600'
                                    }`}>
                                    {obs.type} Impact
                                </span>
                                <p className="text-[11px] font-medium text-darkGrey mt-0.5">{obs.text}</p>
                            </div>
                        </div>
                    ))}
                </section>

                {/* Contradictions / Alerts */}
                {insight.contradictions.length > 0 && (
                    <section className="space-y-3">
                        <h4 className="text-[10px] font-black uppercase text-red-600 tracking-wider">Scientific Discrepancies</h4>
                        {insight.contradictions.map((c, idx) => (
                            <div key={idx} className="p-3 bg-red-50 border border-red-100 rounded-xl flex gap-3">
                                <ShieldAlert size={16} className="text-red-500 shrink-0" />
                                <p className="text-[11px] font-bold text-red-800 leading-tight">
                                    <span className="uppercase">{c.source1} vs {c.source2}:</span> {c.text}
                                </p>
                            </div>
                        ))}
                    </section>
                )}

                {/* Suggested Next Steps */}
                <section className="space-y-3">
                    <h4 className="text-[10px] font-black uppercase text-midGrey tracking-wider">Suggested Actions</h4>
                    {insight.suggested_next_steps.map((step, idx) => (
                        <button key={idx} className="w-full p-4 border border-[#E5E7EB] hover:border-blue-500 text-left rounded-2xl group transition-all hover:bg-blue-50/30">
                            <div className="flex justify-between items-center mb-1">
                                <span className="text-[11px] font-black uppercase text-blue-600">{step.action}</span>
                                <ChevronRight size={14} className="text-gray-300 group-hover:text-blue-600" />
                            </div>
                            <p className="text-[10px] font-medium text-midGrey leading-snug">
                                {step.rationale}
                            </p>
                        </button>
                    ))}
                </section>
            </div>
        </div>
    );
}
