import React, { useState } from 'react';
import { useLabStore } from "../../../store/labStore";

const ATOMS = [
    { symbol: 'C', color: '#2B2B2B', name: 'Carbon' },
    { symbol: 'H', color: '#EDEDED', name: 'Hydrogen' },
    { symbol: 'O', color: '#E53935', name: 'Oxygen' },
    { symbol: 'N', color: '#1E88E5', name: 'Nitrogen' },
    { symbol: 'S', color: '#FBC02D', name: 'Sulfur' },
    { symbol: 'F', color: '#43A047', name: 'Fluorine' },
    { symbol: 'Cl', color: '#2E7D32', name: 'Chlorine' },
    { symbol: 'Br', color: '#8D6E63', name: 'Bromine' },
];

const BONDS = [
    { type: 1, label: 'Single', icon: '—' },
    { type: 2, label: 'Double', icon: '==' },
    { type: 3, label: 'Triple', icon: '≡' },
    { type: 4, label: 'Aromatic', icon: '⏣' },
];

export default function LeftToolDock() {
    const { currentTool, setTool } = useLabStore();
    // In a real implementation, we would store selectedElement/Bond separately in the store
    // For now, we'll just log or use local state for visual feedback
    const [selectedElement, setSelectedElement] = useState('C');
    const [selectedBondData, setSelectedBondData] = useState(1);

    return (
        <div
            className="
                fixed left-4 top-28 bottom-24
                w-16 hover:w-52
                transition-all duration-200 ease-out
                bg-white/90 backdrop-blur
                rounded-2xl shadow-md
                overflow-hidden
                flex flex-col
                group
                z-50
            "
        >
            <div className="flex-1 overflow-y-auto p-2 scrollbar-hide">
                <div className="flex flex-col items-center gap-4 pt-2">
                    {/* Atoms Section */}
                    <div className="w-full">
                        <div className="text-xs font-bold text-gray-400 uppercase mb-2 px-2 opacity-0 group-hover:opacity-100 transition-opacity duration-200 whitespace-nowrap">
                            Atoms
                        </div>
                        <div className="flex flex-col gap-2 items-center group-hover:items-stretch">
                            {ATOMS.map(atom => (
                                <button
                                    key={atom.symbol}
                                    onClick={() => {
                                        setTool('add-atom');
                                        setSelectedElement(atom.symbol);
                                        // Update store appropriately if needed
                                    }}
                                    className={`
                                        flex items-center gap-3 p-1.5 rounded-lg transition-colors
                                        ${selectedElement === atom.symbol && currentTool === 'add-atom' ? 'bg-blue-50 ring-1 ring-blue-500' : 'hover:bg-gray-100'}
                                    `}
                                >
                                    <div
                                        className="w-8 h-8 rounded-full shadow-sm flex items-center justify-center text-xs font-bold shrink-0 border border-black/10"
                                        style={{ backgroundColor: atom.color, color: atom.symbol === 'C' ? 'white' : 'black' }}
                                    >
                                        {atom.symbol}
                                    </div>
                                    <span className="text-sm font-medium text-gray-700 opacity-0 group-hover:opacity-100 transition-opacity duration-150 whitespace-nowrap">
                                        {atom.name}
                                    </span>
                                </button>
                            ))}
                        </div>
                    </div>

                    <div className="w-full h-px bg-gray-200 my-1" />

                    {/* Bonds Section */}
                    <div className="w-full mb-4">
                        <div className="text-xs font-bold text-gray-400 uppercase mb-2 px-2 opacity-0 group-hover:opacity-100 transition-opacity duration-200 whitespace-nowrap">
                            Bonds
                        </div>
                        <div className="flex flex-col gap-2 items-center group-hover:items-stretch">
                            {BONDS.map(bond => (
                                <button
                                    key={bond.type}
                                    onClick={() => {
                                        setTool('add-bond');
                                        setSelectedBondData(bond.type);
                                    }}
                                    className={`
                                        flex items-center gap-3 p-1.5 rounded-lg transition-colors
                                        ${selectedBondData === bond.type && currentTool === 'add-bond' ? 'bg-blue-50 ring-1 ring-blue-500' : 'hover:bg-gray-100'}
                                    `}
                                >
                                    <div className="w-8 h-8 rounded-md bg-white border border-gray-200 flex items-center justify-center shrink-0 font-mono text-xs">
                                        {bond.icon}
                                    </div>
                                    <span className="text-sm font-medium text-gray-700 opacity-0 group-hover:opacity-100 transition-opacity duration-150 whitespace-nowrap">
                                        {bond.label}
                                    </span>
                                </button>
                            ))}
                        </div>
                    </div>
                </div>
            </div>

            {/* Template Launcher / Footer */}
            <div className="p-2 border-t border-gray-100 bg-gray-50/50">
                <button className="w-full flex items-center gap-3 p-1.5 rounded-lg hover:bg-gray-100 text-gray-600">
                    <div className="w-8 h-8 flex items-center justify-center shrink-0">
                        📁
                    </div>
                    <span className="text-sm font-medium opacity-0 group-hover:opacity-100 transition-opacity duration-150 whitespace-nowrap">
                        Templates
                    </span>
                </button>
            </div>
        </div>
    );
}
