import React from 'react';
import { useLabStore } from "../../../store/labStore";

export default function RightInspector() {
    // In a real implementation, we would check selectedAtomId/selectedBondId from store
    // For skeleton, we'll assume it's hidden unless explicitly shown for demo
    // const { selectedEntity } = useLabStore();
    const selectedEntity = null; // Mock for now

    // Check if we show it
    const isVisible = !!selectedEntity;

    return (
        <div
            className={`
                fixed right-4 top-28 bottom-24
                w-80
                bg-white/95 backdrop-blur
                rounded-2xl shadow-lg
                transition-transform duration-300 ease-in-out
                flex flex-col overflow-hidden
                z-50
                ${isVisible ? 'translate-x-0 opacity-100' : 'translate-x-[110%] opacity-0 pointer-events-none'}
            `}
        >
            <div className="p-4 border-b border-gray-100 flex justify-between items-center bg-gray-50/50">
                <h3 className="font-semibold text-gray-800">Inspector</h3>
                <span className="text-xs text-gray-400 uppercase tracking-wider">Properties</span>
            </div>

            <div className="flex-1 p-4 overflow-y-auto">
                <div className="flex items-center justify-center h-full text-gray-400 text-sm italic">
                    Select an item to view properties
                </div>
            </div>
        </div>
    );
}
