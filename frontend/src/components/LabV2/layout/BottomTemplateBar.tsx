
import React from 'react';
import { useLabStore } from "../../../store/labStore";
import { BENZENE, CYCLOHEXANE, CYCLOPROPANE, CHAIN_4 } from "../../../utils/defaultMolecules";
import { Hexagon, Triangle, Circle, Share2 } from "lucide-react";

const TEMPLATES = [
    { id: 'benzene', name: 'Benzene', data: BENZENE, icon: Hexagon },
    { id: 'cyclohexane', name: 'Cyclohexane', data: CYCLOHEXANE, icon: Circle },
    { id: 'chain', name: 'Chain', data: CHAIN_4, icon: Share2 },
    { id: 'triangle', name: 'Cyclopropane', data: CYCLOPROPANE, icon: Triangle },
];

export default function BottomTemplateBar() {
    const { loadMolecule } = useLabStore();

    return (
        <div className="w-full h-full bg-white border-t border-gray-100 px-6 flex items-center gap-6 overflow-x-auto">
            <span className="text-[10px] uppercase font-bold text-gray-400 tracking-wider shrink-0">Templates</span>

            <div className="h-8 w-px bg-gray-200 shrink-0" />

            <div className="flex items-center gap-4">
                {TEMPLATES.map(t => {
                    const Icon = t.icon;
                    return (
                        <button
                            key={t.id}
                            onClick={() => {
                                if (t.data) {
                                    // Map legacy template data to new Core structure
                                    const d = t.data as any; // Legacy structure
                                    const atomMap = new Map<string, string>(); // Legacy ID -> New ID

                                    const atoms = d.atoms.map((a: any) => {
                                        const newId = `atom-${Math.random().toString(36).substr(2, 9)}`;
                                        atomMap.set(a.id, newId);

                                        // Handle both old schema (x,y,z) and new schema (position array)
                                        const pos = Array.isArray(a.position)
                                            ? { x: a.position[0], y: a.position[1], z: a.position[2] }
                                            : { x: a.x, y: a.y, z: a.z || 0 };

                                        return {
                                            id: newId,
                                            element: a.element,
                                            position: pos, // Map to object
                                            charge: 0
                                        };
                                    });

                                    const bonds = d.bonds.map((b: any) => ({
                                        id: `bond-${Math.random().toString(36).substr(2, 9)}`,
                                        from: atomMap.get(b.from || b.a)!,
                                        to: atomMap.get(b.to || b.b)!,
                                        order: b.order || 1
                                    }));

                                    loadMolecule({
                                        id: `mol-${Date.now()}`,
                                        atoms,
                                        bonds,
                                        metadata: {}
                                    });
                                }
                            }}
                            className="flex items-center gap-2 px-3 py-1.5 rounded-lg hover:bg-gray-50 text-gray-600 transition-colors border border-transparent hover:border-gray-200 group"
                        >
                            <Icon size={16} className="text-gray-400 group-hover:text-blue-500 transition-colors" />
                            <span className="text-sm font-medium">{t.name}</span>
                        </button>
                    )
                })}
            </div>
        </div>
    );
}
