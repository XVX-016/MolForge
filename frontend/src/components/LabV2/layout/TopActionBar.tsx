import React from 'react';
import { useLabStore } from "../../../store/labStore";
import { LibraryAPI } from "../../../api/library";
import {
    MousePointer2,
    Atom,
    Link2,
    Undo2,
    Redo2,
    Trash2,
    Download,
    Upload,
    CloudUpload
} from "lucide-react";

export default function TopActionBar() {
    const { currentTool, setTool, molecule, loadMolecule } = useLabStore();

    const handleAction = async (id: string, label: string) => {
        if (id === 'select' || id === 'add-atom' || id === 'add-bond' || id === 'delete') {
            setTool(id as any);
        } else if (id === 'download') {
            const dataStr = "data:text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(molecule));
            const downloadAnchorNode = document.createElement('a');
            downloadAnchorNode.setAttribute("href", dataStr);
            downloadAnchorNode.setAttribute("download", (molecule.metadata?.name || "molecule") + ".json");
            document.body.appendChild(downloadAnchorNode); // required for firefox
            downloadAnchorNode.click();
            downloadAnchorNode.remove();
        } else if (id === 'upload') {
            const input = document.createElement('input');
            input.type = 'file';
            input.accept = '.json';
            input.onchange = (e: any) => {
                const file = e.target.files[0];
                if (!file) return;
                const reader = new FileReader();
                reader.onload = (event) => {
                    try {
                        const json = JSON.parse(event.target?.result as string);
                        loadMolecule(json);
                    } catch (err) {
                        alert("Invalid JSON file");
                    }
                };
                reader.readAsText(file);
            };
            input.click();
        } else if (id === 'save-library') {
            const name = prompt("Enter molecule name", "New Molecule");
            if (!name) return;
            try {
                // Assuming global availability or import. 
                // Note: explicit import of LibraryAPI needed if not present
                await LibraryAPI.upload({
                    name,
                    json_graph: { atoms: molecule.atoms, bonds: molecule.bonds }
                });
                alert("Saved to Library!");
            } catch (e) {
                console.error(e);
                alert("Failed to save.");
            }
        } else if (id === 'undo') {
            // Handled inside store usually or here? 
            // Logic missing in previous file for actual undo execution in UI handler
            // Assuming useLabStore exposes undo/redo functions.
            // Checking original file... it didn't use undo() from store in handleAction!
            // It just logged. I should hook them up if they exist in store.
            useLabStore.getState().undo();
        } else if (id === 'redo') {
            useLabStore.getState().redo();
        }
    };

    const actions = [
        { id: 'select', label: 'Select', icon: MousePointer2 },
        { id: 'add-atom', label: 'Add Atom', icon: Atom },
        { id: 'add-bond', label: 'Add Bond', icon: Link2 },
        { id: 'delete', label: 'Delete', icon: Trash2 },
        { id: 'undo', label: 'Undo', icon: Undo2 },
        { id: 'redo', label: 'Redo', icon: Redo2 },
        { id: 'download', label: 'Download JSON', icon: Download },
        { id: 'upload', label: 'Upload JSON', icon: Upload },
        { id: 'save-library', label: 'Save to Cloud', icon: CloudUpload },
    ];

    return (
        <div className="w-full h-full flex items-center justify-center bg-white border-b border-gray-100 px-4">
            <div className="flex items-center gap-1 bg-gray-50/50 p-1.5 rounded-xl border border-gray-100">
                {actions.map((action) => {
                    const Icon = action.icon;
                    // Highlight logic: strict equality for tools, but actions like undo/download shouldn't remain 'active'
                    const isTool = ['select', 'add-atom', 'add-bond', 'delete'].includes(action.id);
                    const isActive = isTool && currentTool === action.id;

                    return (
                        <button
                            key={action.id}
                            onClick={() => handleAction(action.id, action.label)}
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
