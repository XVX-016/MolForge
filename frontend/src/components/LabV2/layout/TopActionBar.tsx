import React from 'react';
import { useLabStore } from "../../../store/labStore";

export default function TopActionBar() {
    const { currentTool, setTool } = useLabStore();

    const actions = [
        { id: 'select', label: 'Select', icon: '🖱️' },
        { id: 'add-atom', label: 'Add Atom', icon: '⚛️' },
        { id: 'add-bond', label: 'Add Bond', icon: '🔗' },
        { id: 'autobond', label: 'AutoBond', icon: '✨' },
        { id: 'optimize', label: 'Optimize', icon: '⚡' },
        { id: 'undo', label: 'Undo', icon: '↩️' },
        { id: 'redo', label: 'Redo', icon: '↪️' },
        { id: 'save', label: 'Save', icon: '💾' },
    ];

    return (
        <div className="fixed top-20 left-1/2 -translate-x-1/2 flex gap-2 px-4 py-2 bg-white/80 backdrop-blur-md rounded-full shadow-lg z-50">
            {actions.map((action) => (
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
                        p-2 rounded-full transition-all duration-200
                        ${currentTool === action.id ? 'bg-black text-white shadow-md scale-105' : 'text-gray-600 hover:bg-gray-100'}
                    `}
                    title={action.label}
                >
                    <span className="text-lg">{action.icon}</span>
                </button>
            ))}
        </div>
    );
}
