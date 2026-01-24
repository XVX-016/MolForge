import { useStudioStore } from '../../store/studioStore';
import { motion } from 'framer-motion';
import { Database, Check, History, Zap } from 'lucide-react';
import { RadarChart } from './RadarChart';

export default function PropertiesPanel() {
    const {
        properties,
        isComputingProperties,
        commitDraft,
        isSyncing,
        activeVersionId,
        activeVersionIndex,
        isDirty,
        canCommit,
        analysis
    } = useStudioStore();

    if (!properties) return null;

    const items = [
        { label: 'MW', value: properties.molecular_weight, unit: 'g/mol' },
        { label: 'LogP', value: properties.logp, unit: '' },
        { label: 'TPSA', value: properties.tpsa, unit: 'Å²' },
        { label: 'HBD', value: properties.hbd, unit: '' },
        { label: 'HBA', value: properties.hba, unit: '' },
        { label: 'RotB', value: properties.rotatable_bonds, unit: '' },
        { label: 'QED', value: properties.qed, unit: '', highlight: true },
        { label: 'Cmplx', value: properties.complexity, unit: '' },
    ];

    return (
        <div className="bg-white/90 backdrop-blur-md border border-gray-200 rounded-xl p-4 shadow-sm w-full font-sans">
            <div className="flex items-center justify-between mb-4">
                <div className="flex items-center gap-2">
                    <History className="w-3 h-3 text-gray-400" />
                    <h3 className="text-[10px] font-black uppercase tracking-[0.2em] text-gray-400">
                        {activeVersionId ? `Version ${activeVersionIndex}` : 'Draft Snapshot'}
                    </h3>
                </div>
                {isComputingProperties && (
                    <motion.div
                        animate={{ rotate: 360 }}
                        transition={{ repeat: Infinity, duration: 1, ease: "linear" }}
                        className="w-3 h-3 border border-blue-500 border-t-transparent rounded-full"
                    />
                )}
            </div>

            <div className="grid grid-cols-2 gap-x-6 gap-y-3">
                {items.map((item) => (
                    <div key={item.label} className="flex flex-col">
                        <span className={`text-[9px] font-bold uppercase tracking-wider ${item.highlight ? 'text-blue-500' : 'text-gray-400'}`}>
                            {item.label}
                        </span>
                        <div className="flex items-baseline gap-1">
                            <span className={`text-sm font-semibold ${item.highlight ? 'text-blue-600' : 'text-gray-900'}`}>
                                {item.value}
                            </span>
                            {item.unit && <span className="text-[10px] text-gray-400">{item.unit}</span>}
                        </div>
                    </div>
                ))}
            </div>

            {/* Radar Analysis (Phase 5B) */}
            {analysis?.radar && (
                <div className="mt-4 pt-4 border-t border-gray-100">
                    <div className="flex items-center gap-2 mb-2">
                        <Zap className="w-3 h-3 text-amber-500" />
                        <span className="text-[10px] font-black uppercase tracking-wider text-gray-500">Molecular Profile</span>
                    </div>
                    <div className="bg-gray-50/50 rounded-xl p-2 border border-blue-50/50">
                        <RadarChart scores={analysis.radar} />
                    </div>
                </div>
            )}

            {/* QED Visualization (Phase 5B) */}
            <div className="mt-4 px-1">
                <div className="flex items-center justify-between mb-1">
                    <span className="text-[8px] font-black uppercase text-gray-400">Drug-likeness (QED)</span>
                    <span className="text-[10px] font-bold text-gray-700">{properties.qed}</span>
                </div>
                <div className="h-1.5 w-full bg-gray-100 rounded-full overflow-hidden">
                    <motion.div
                        initial={{ width: 0 }}
                        animate={{ width: `${properties.qed * 100}%` }}
                        className={`h-full rounded-full ${properties.qed > 0.6 ? 'bg-emerald-500' :
                            properties.qed > 0.4 ? 'bg-amber-400' : 'bg-rose-500'
                            }`}
                    />
                </div>
            </div>

            {/* Identifiers (Auditability) */}
            <div className="mt-4 space-y-2">
                <div className="flex flex-col gap-0.5">
                    <span className="text-[8px] font-bold text-gray-400 uppercase tracking-wider">InChIKey</span>
                    <span className="text-[9px] font-mono break-all text-gray-600 bg-gray-50 px-1.5 py-0.5 rounded border border-gray-100">
                        {properties.inchikey}
                    </span>
                </div>
                <div className="flex flex-col gap-0.5">
                    <span className="text-[8px] font-bold text-gray-400 uppercase tracking-wider">Canonical SMILES</span>
                    <span className="text-[9px] font-mono break-all text-gray-500 line-clamp-2">
                        {properties.canonical_smiles}
                    </span>
                </div>
            </div>

            <div className="mt-4 pt-4 border-t border-gray-100 flex items-center justify-between">
                <div className="flex flex-col">
                    <span className="text-[9px] font-bold text-gray-400 uppercase tracking-wider">Lipinski</span>
                    <span className={`text-sm font-semibold ${properties.lipinski_violations > 0 ? 'text-amber-500' : 'text-emerald-500'}`}>
                        {properties.lipinski_violations > 0 ? `${properties.lipinski_violations} Violations` : 'Passed'}
                    </span>
                </div>
                <div className="flex flex-col text-right">
                    <span className="text-[9px] font-bold text-gray-400 uppercase tracking-wider">Formula</span>
                    <span className="text-xs font-mono font-bold text-gray-700">{properties.formula}</span>
                </div>
            </div>

            {/* Commit Button */}
            <button
                onClick={() => commitDraft()}
                disabled={isSyncing || !canCommit || !isDirty}
                className={`mt-4 w-full flex items-center justify-center gap-2 py-2 rounded-lg text-xs font-bold transition-all ${!isDirty && activeVersionId
                    ? 'bg-emerald-50 text-emerald-600 border border-emerald-100'
                    : !canCommit || !isDirty
                        ? 'bg-gray-100 text-gray-400 cursor-not-allowed border border-gray-200'
                        : 'bg-blue-600 text-white hover:bg-blue-700 shadow-sm'
                    }`}
            >
                {isSyncing ? (
                    <motion.div
                        animate={{ rotate: 360 }}
                        transition={{ repeat: Infinity, duration: 1, ease: "linear" }}
                        className="w-3 h-3 border-2 border-white/30 border-t-white rounded-full"
                    />
                ) : !isDirty && activeVersionId ? (
                    <><Check className="w-3 h-3 text-emerald-500" /> Versioned Snapshot</>
                ) : (
                    <><Database className="w-3 h-3" /> Commit Snapshot</>
                )}
            </button>

            <p className="mt-3 text-[8px] text-gray-400 italic text-center">
                Molecular Intelligence locked to RDKit Identity
            </p>
        </div>
    );
}
