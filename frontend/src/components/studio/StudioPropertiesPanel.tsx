
import { useStudioStore } from '../../store/studioStore';
import { AlertCircle, Zap, Activity, Info, TrendingDown, TrendingUp } from 'lucide-react';
import { RadarChart } from './RadarChart';

export default function StudioPropertiesPanel() {
    const { status, dashboard } = useStudioStore();

    if (!dashboard) return null;

    const { baseline, proposal, alerts, property_delta, radar } = dashboard;

    const properties = proposal ? proposal.properties : baseline.properties;

    const propertyConfigs = [
        { key: 'molecular_weight', label: 'MW', unit: 'g/mol' },
        { key: 'logp', label: 'LogP', unit: '' },
        { key: 'tpsa', label: 'TPSA', unit: 'Å²' },
        { key: 'qed', label: 'QED', unit: '' },
    ];

    return (
        <div className="flex flex-col p-6 space-y-8">
            {/* 1. Alerts & Audits */}
            <section>
                <h3 className="text-[10px] font-black text-gray-400 uppercase tracking-widest mb-4 flex items-center gap-2">
                    <AlertCircle size={12} className="text-red-500" /> Clinical Audit
                </h3>
                <div className="space-y-2">
                    {alerts.map((alert, i) => (
                        <div key={i} className={`p-3 rounded-lg border flex gap-3 ${alert.severity === 'high' ? 'bg-red-50 border-red-100' : 'bg-amber-50 border-amber-100'
                            }`}>
                            <AlertCircle size={14} className={alert.severity === 'high' ? 'text-red-600' : 'text-amber-600'} />
                            <div className="flex flex-col">
                                <span className={`text-[10px] font-bold uppercase tracking-tight ${alert.severity === 'high' ? 'text-red-700' : 'text-amber-700'
                                    }`}>{alert.id}</span>
                                <p className="text-[10px] text-gray-600 leading-tight">{alert.description}</p>
                            </div>
                        </div>
                    ))}
                    {alerts.length === 0 && (
                        <div className="py-2 text-[11px] text-green-600 flex items-center gap-2">
                            <Zap size={12} /> No structural alerts detected.
                        </div>
                    )}
                </div>
            </section>

            {/* 2. RDKit Properties & Deltas */}
            <section>
                <h3 className="text-[10px] font-black text-gray-400 uppercase tracking-widest mb-4 flex items-center gap-2">
                    <Activity size={12} className="text-blue-600" /> RDKit Analytics
                </h3>
                <div className="grid grid-cols-2 gap-4">
                    {propertyConfigs.map((cfg) => {
                        const delta = property_delta[cfg.key];
                        const val = properties[cfg.key];
                        return (
                            <div key={cfg.key} className="p-4 rounded-xl border border-[#E5E7EB] bg-white shadow-sm relative overflow-hidden">
                                <span className="text-[9px] font-black text-gray-400 uppercase tracking-wider">{cfg.label}</span>
                                <div className="flex items-baseline gap-1 mt-1">
                                    <span className="text-lg font-bold text-gray-900">{val}</span>
                                    <span className="text-[9px] text-gray-400 font-medium">{cfg.unit}</span>
                                </div>
                                {delta !== undefined && delta !== 0 && (
                                    <div className={`mt-2 flex items-center gap-1 text-[9px] font-bold ${delta < 0 ? 'text-green-600' : 'text-red-500'
                                        }`}>
                                        {delta < 0 ? <TrendingDown size={10} /> : <TrendingUp size={10} />}
                                        {delta > 0 ? '+' : ''}{delta}
                                    </div>
                                )}
                            </div>
                        );
                    })}
                </div>
            </section>

            {/* 3. Radar Comparison */}
            <section>
                <div className="flex items-center justify-between mb-4">
                    <h3 className="text-[10px] font-black text-gray-400 uppercase tracking-widest flex items-center gap-2">
                        <Activity size={12} className="text-purple-600" /> Multi-Objective Profile
                    </h3>
                    <div className="flex gap-3">
                        <div className="flex items-center gap-1">
                            <div className="w-1.5 h-1.5 rounded-full bg-gray-300" />
                            <span className="text-[8px] font-bold text-gray-400 uppercase">Baseline</span>
                        </div>
                        <div className="flex items-center gap-1">
                            <div className="w-1.5 h-1.5 rounded-full bg-blue-500" />
                            <span className="text-[8px] font-bold text-blue-500 uppercase">Proposal</span>
                        </div>
                    </div>
                </div>
                <div className="bg-[#F9FAFB] border border-[#E5E7EB] rounded-2xl p-4 aspect-square">
                    <RadarChart baseline={radar.baseline} proposal={radar.proposal} />
                </div>
            </section>

            {/* 4. Persistence Controls (Fixed Bottom) */}
            <div className="mt-auto pt-6 border-t border-[#E5E7EB] space-y-3">
                {status === 'PROPOSED' && (
                    <div className="flex gap-2">
                        <button
                            className="flex-1 bg-blue-600 text-white py-3 rounded-xl text-xs font-black uppercase tracking-widest hover:bg-blue-700 transition-all shadow-lg shadow-blue-200"
                            onClick={() => { }} // acceptProposal
                        >
                            Commit Proposal
                        </button>
                        <button
                            className="px-4 py-3 border border-[#E5E7EB] text-gray-400 rounded-xl hover:bg-red-50 hover:text-red-500 hover:border-red-100 transition-all"
                            onClick={() => { }} // rejectProposal
                        >
                            ✕
                        </button>
                    </div>
                )}
                <div className="flex items-center justify-center gap-2 opacity-40">
                    <Info size={10} />
                    <span className="text-[8px] font-bold uppercase tracking-tighter">Identity-First Persistence Active</span>
                </div>
            </div>
        </div>
    );
}
