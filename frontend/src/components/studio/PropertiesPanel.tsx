import { useStudioStore } from '../../store/studioStore';
import { motion } from 'framer-motion';

export default function PropertiesPanel() {
    const { properties, isComputingProperties } = useStudioStore();

    if (!properties) return null;

    const items = [
        { label: 'MW', value: properties.molecular_weight, unit: 'g/mol' },
        { label: 'LogP', value: properties.logp, unit: '' },
        { label: 'TPSA', value: properties.tpsa, unit: 'Å²' },
        { label: 'HBD', value: properties.hbd, unit: '' },
        { label: 'HBA', value: properties.hba, unit: '' },
        { label: 'RotB', value: properties.rotatable_bonds, unit: '' },
    ];

    return (
        <div className="bg-white/90 backdrop-blur-md border border-gray-200 rounded-xl p-4 shadow-sm w-full">
            <div className="flex items-center justify-between mb-3">
                <h3 className="text-[10px] font-black uppercase tracking-[0.2em] text-gray-400">
                    Molecular Analysis
                </h3>
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
                        <span className="text-[9px] font-bold text-gray-400 uppercase tracking-wider">{item.label}</span>
                        <div className="flex items-baseline gap-1">
                            <span className="text-sm font-semibold text-gray-900">{item.value}</span>
                            {item.unit && <span className="text-[10px] text-gray-400">{item.unit}</span>}
                        </div>
                    </div>
                ))}
            </div>

            <div className="mt-4 pt-3 border-t border-gray-100 flex items-center justify-between">
                <div className="flex flex-col">
                    <span className="text-[9px] font-bold text-gray-400 uppercase tracking-wider">Lipinski Violations</span>
                    <span className={`text-sm font-semibold ${properties.lipinski_violations > 0 ? 'text-amber-500' : 'text-emerald-500'}`}>
                        {properties.lipinski_violations}
                    </span>
                </div>
                <div className="flex flex-col text-right">
                    <span className="text-[9px] font-bold text-gray-400 uppercase tracking-wider">Formula</span>
                    <span className="text-xs font-mono font-bold text-gray-700">{properties.formula}</span>
                </div>
            </div>

            <p className="mt-3 text-[8px] text-gray-400 italic text-center">
                Deterministic RDKit Computation
            </p>
        </div>
    );
}
