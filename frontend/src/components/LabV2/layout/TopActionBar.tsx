import React from 'react';
import { useLabStore } from "../../../store/labStore";
import {
    MousePointer2,
    Atom,
    Link2,
    Sparkles,
    Zap,
    Undo2,
    Redo2,
    Save
} from "lucide-react";

export default function TopActionBar() {
    const { currentTool, setTool } = useLabStore();

    const actions = [
        { id: 'select', label: 'Select', icon: MousePointer2 },
        { id: 'add-atom', label: 'Add Atom', icon: Atom },
        { id: 'add-bond', label: 'Add Bond', icon: Link2 },
        { id: 'autobond', label: 'AutoBond', icon: Sparkles },
        { id: 'optimize', label: 'Optimize', icon: Zap },
        { id: 'undo', label: 'Undo', icon: Undo2 },
        { id: 'redo', label: 'Redo', icon: Redo2 },
        { id: 'save', label: 'Save', icon: Save },
    ];

    return (
        <div className="w-full h-full flex items-center justify-center bg-white border-b border-gray-100 px-4">
            <div className="flex items-center gap-1 bg-gray-50/50 p-1.5 rounded-xl border border-gray-100">
                {actions.map((action) => {
                    const Icon = action.icon;
                    const isActive = currentTool === action.id;

                    return (
                        <button
                            key={action.id}
                            onClick={() => {
                                if (action.id === 'select' || action.id === 'add-atom' || action.id === 'add-bond') {
                                    setTool(action.id as any);
                                } else {
                                    console.log(`${action.label} triggered`);
                                }
                            }}
                            className={`
                                w-9 h-9 rounded-lg flex items-center justify-center transition-all duration-200
                                ${isActive
                                    ? 'bg-white text-blue-600 shadow-sm ring-1 ring-gray-200'
                                    : 'text-gray-500 hover:bg-gray-100 hover:text-gray-900'
                                }
                            `}
                            title={action.label}
                        >
                            <Icon size={18} strokeWidth={2} />
                        </button>
                    );
                })}
            </div>
        </div>
    );
}
