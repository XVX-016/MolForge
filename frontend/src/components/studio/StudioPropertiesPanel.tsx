
import { useStudioStore } from '../../store/studioStore';
import { AlertCircle, Activity, Info, TrendingDown, TrendingUp, CheckCircle, XCircle, ShieldCheck } from 'lucide-react';
import { RadarChart } from './RadarChart';

export default function AuditPanel() {
    const { status, dashboard, acceptProposal, rejectProposal } = useStudioStore();

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
        <div className="flex flex-col h-full bg-white border-l border-[#E5E7EB]">
            <div className="p-6 overflow-y-auto flex-1 space-y-8">
                {/* 1. Alerts & Audits */}
                <section>
                    <div className="flex items-center justify-between mb-4">
                        <h3 className="text-[10px] font-black text-gray-400 uppercase tracking-widest flex items-center gap-2">
                            <AlertCircle size={12} className="text-red-500" /> Clinical Audit
                        </h3>
                        {alerts.length > 0 && (
                            <span className="text-[8px] font-black px-1.5 py-0.5 bg-red-50 text-red-600 rounded uppercase tracking-tighter">
                                {alerts.filter(a => a.severity === 'high').length} High Severity
                            </span>
                        )}
                    </div>

                    <div className="space-y-2">
                        {alerts.map((alert, i) => (
                            <div key={i} className={`p-3 rounded-lg border flex gap-3 transition-all ${alert.severity === 'high'
                                ? 'bg-red-50/50 border-red-100 hover:bg-red-50'
                                : 'bg-amber-50/50 border-amber-100 hover:bg-amber-50'
                                }`}>
                                <AlertCircle size={14} className={alert.severity === 'high' ? 'text-red-500 mt-0.5' : 'text-amber-500 mt-0.5'} />
                                <div className="flex flex-col">
                                    <span className={`text-[9px] font-black uppercase tracking-tight ${alert.severity === 'high' ? 'text-red-700' : 'text-amber-700'
                                        }`}>{alert.id}</span>
                                    <p className="text-[10px] text-gray-600 leading-tight mt-0.5">{alert.description}</p>
                                </div>
                            </div>
                        ))}
                        {alerts.length === 0 && (
                            <div className="py-4 rounded-xl border border-dashed border-green-200 bg-green-50/20 flex flex-col items-center justify-center text-center">
                                <CheckCircle size={16} className="text-green-500 mb-2 opacity-50" />
                                <p className="text-[10px] text-green-700 font-medium tracking-tight">Structurally Clean</p>
                                <p className="text-[8px] text-green-600/60 uppercase font-black mt-1">No PAINS/Lipinski Violations</p>
                            </div>
                        )}
                    </div>
                </section>

                {/* 2. RDKit Analytics */}
                <section>
                    <h3 className="text-[10px] font-black text-gray-400 uppercase tracking-widest mb-4 flex items-center gap-2">
                        <Activity size={12} className="text-blue-600" /> RDKit Analytics
                    </h3>
                    <div className="grid grid-cols-2 gap-3">
                        {propertyConfigs.map((cfg) => {
                            const delta = property_delta[cfg.key];
                            const val = properties[cfg.key];
                            return (
                                <div key={cfg.key} className="p-3 rounded-xl border border-[#F3F4F6] bg-white shadow-[0_1px_2px_rgba(0,0,0,0.02)] transition-all hover:border-blue-100 group">
                                    <div className="flex justify-between items-start mb-1">
                                        <span className="text-[8px] font-black text-gray-400 uppercase tracking-widest group-hover:text-blue-400 transition-colors">{cfg.label}</span>
                                        {delta !== undefined && delta !== 0 && (
                                            <div className={`flex items-center gap-0.5 text-[8px] font-black ${delta < 0 ? 'text-green-600' : 'text-red-500'
                                                }`}>
                                                {delta < 0 ? <TrendingDown size={8} /> : <TrendingUp size={8} />}
                                                {Math.abs(delta)}
                                            </div>
                                        )}
                                    </div>
                                    <div className="flex items-baseline gap-1">
                                        <span className="text-sm font-bold text-gray-900">{typeof val === 'number' ? val.toFixed(2) : val}</span>
                                        <span className="text-[8px] text-gray-300 font-bold uppercase">{cfg.unit}</span>
                                    </div>
                                </div>
                            );
                        })}
                    </div>
                </section>

                {/* 3. Multi-Objective Profile */}
                <section>
                    <div className="flex items-center justify-between mb-4">
                        <h3 className="text-[10px] font-black text-gray-400 uppercase tracking-widest flex items-center gap-2">
                            <Activity size={12} className="text-purple-600" /> Objective Map
                        </h3>
                        <div className="flex gap-2">
                            <div className="flex items-center gap-1.5 px-2 py-0.5 rounded-full bg-slate-50 border border-slate-100">
                                <div className="w-1.5 h-1.5 rounded-full bg-gray-300" />
                                <span className="text-[7px] font-black text-gray-400 uppercase tracking-tighter">v.Baseline</span>
                            </div>
                            <div className="flex items-center gap-1.5 px-2 py-0.5 rounded-full bg-blue-50 border border-blue-100">
                                <div className="w-1.5 h-1.5 rounded-full bg-blue-500" />
                                <span className="text-[7px] font-black text-blue-500 uppercase tracking-tighter">v.Proposal</span>
                            </div>
                        </div>
                    </div>
                    <div className="bg-slate-50/50 border border-[#F3F4F6] rounded-2xl p-4 aspect-square flex items-center justify-center">
                        <RadarChart baseline={radar.baseline} proposal={radar.proposal} />
                    </div>
                </section>
            </div>

            {/* 4. Persistence & Sign-off */}
            <div className="p-6 border-t border-[#E5E7EB] bg-white space-y-4 shadow-[0_-4px_12px_rgba(0,0,0,0.02)]">
                {status === 'PROPOSED' ? (
                    <div className="space-y-3">
                        <div className="bg-blue-50/50 border border-blue-100 rounded-xl p-3">
                            <div className="flex items-center gap-2 mb-1">
                                <Info size={12} className="text-blue-600" />
                                <span className="text-[9px] font-bold text-blue-700 uppercase">Awaiting Chemist Sign-off</span>
                            </div>
                            <p className="text-[9px] text-blue-600 leading-tight italic">
                                Verification of atomic changes and property deltas required for commit.
                            </p>
                        </div>
                        <div className="flex gap-2">
                            <button
                                className="flex-1 bg-blue-600 text-white py-3 rounded-xl text-[10px] font-black uppercase tracking-widest hover:bg-blue-700 transition-all shadow-lg shadow-blue-200 active:scale-95 flex items-center justify-center gap-2"
                                onClick={acceptProposal}
                            >
                                <CheckCircle size={14} /> Commit Proposal
                            </button>
                            <button
                                className="px-4 py-3 border border-[#E5E7EB] text-gray-400 rounded-xl hover:bg-red-50 hover:text-red-500 hover:border-red-100 transition-all active:scale-95"
                                onClick={rejectProposal}
                            >
                                <XCircle size={14} />
                            </button>
                        </div>
                    </div>
                ) : status === 'COMMITTED' ? (
                    <div className="py-8 flex flex-col items-center justify-center animate-in fade-in zoom-in duration-300">
                        <div className="w-12 h-12 rounded-full bg-green-50 flex items-center justify-center mb-3">
                            <CheckCircle size={24} className="text-green-500" />
                        </div>
                        <span className="text-[11px] font-bold text-gray-900 uppercase tracking-widest">Version Committed</span>
                        <span className="text-[9px] text-gray-500 font-mono mt-1 uppercase">Audit Log Anchored</span>
                    </div>
                ) : (
                    <div className="flex flex-col items-center justify-center py-4 bg-slate-50/50 rounded-xl border border-dashed border-[#E5E7EB] opacity-60">
                        <span className="text-[9px] font-bold text-gray-400 uppercase tracking-widest">Baseline Monitoring</span>
                    </div>
                )}

                <div className="flex items-center justify-center gap-2">
                    <ShieldCheck size={10} className="text-gray-300" />
                    <span className="text-[8px] font-black text-gray-300 uppercase tracking-widest">21 CFR Part 11 Compliant Audit Log</span>
                </div>
            </div>
        </div>
    );
}
