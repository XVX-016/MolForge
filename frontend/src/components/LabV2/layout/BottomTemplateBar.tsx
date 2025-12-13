import React from 'react';
import { useLabStore } from "../../../store/labStore";
import { BENZENE } from "../../../utils/defaultMolecules";

const TEMPLATES = [
    { id: 'benzene', name: 'Benzene', data: BENZENE, icon: '⏣' },
    { id: 'cyclohexane', name: 'Cyclohexane', data: null, icon: '⬡' }, // Mock data
    { id: 'chain', name: 'Chain', data: null, icon: '〰️' },
    { id: 'group', name: 'Func. Group', data: null, icon: '🧪' },
];

export default function BottomTemplateBar() {
    const { loadMolecule } = useLabStore();

    return (
        <div className="fixed bottom-6 left-1/2 -translate-x-1/2 flex gap-4 px-6 py-3 bg-white/90 backdrop-blur rounded-2xl shadow-md z-50 transition-transform hover:scale-105">
            {TEMPLATES.map(t => (
                <button
                    key={t.id}
                    onClick={() => {
                        console.log(`Add ${t.name}`);
                        // Logic to add template would go here
                    }}
                    className="flex flex-col items-center gap-1 group min-w-[60px]"
                >
                    <div className="w-10 h-10 rounded-lg bg-gray-50 border border-gray-200 flex items-center justify-center text-lg group-hover:bg-blue-50 group-hover:border-blue-200 transition-colors shadow-sm">
                        {t.icon}
                    </div>
                    <span className="text-[10px] uppercase font-bold text-gray-500 group-hover:text-blue-600 tracking-wide">
                        {t.name}
                    </span>
                </button>
            ))}
        </div>
    );
}
