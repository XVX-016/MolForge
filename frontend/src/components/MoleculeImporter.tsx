import React, { useRef } from 'react';
import * as OCL from 'openchemlib';
import { nanoid } from 'nanoid';
import { Upload } from 'lucide-react';
import { supabase } from '../supabase';

interface MoleculeImporterProps {
    userId: string | null;
    onImportSuccess: () => void;
}

export default function MoleculeImporter({ userId, onImportSuccess }: MoleculeImporterProps) {
    const fileInputRef = useRef<HTMLInputElement>(null);

    const handleFileChange = async (e: React.ChangeEvent<HTMLInputElement>) => {
        const file = e.target.files?.[0];
        if (!file || !userId) return;

        try {
            const text = await file.text();
            let mol: any; // OCL.Molecule

            if (file.name.endsWith('.mol') || file.name.endsWith('.sdf')) {
                mol = OCL.Molecule.fromMolfile(text);
            } else if (file.name.endsWith('.smi') || file.name.endsWith('.smiles') || file.name.endsWith('.txt')) {
                mol = OCL.Molecule.fromSmiles(text.trim());
            } else {
                alert('Unsupported file format. Use .mol, .sdf, or .smiles');
                return;
            }
            const atomCount = mol.getAllAtoms();
            const bondCount = mol.getAllBonds();

            const atoms = [];
            const bonds = [];
            const atomIdMap = new Map<number, string>(); // Index -> ID

            for (let i = 0; i < atomCount; i++) {
                const id = nanoid();
                atomIdMap.set(i, id);

                // OCL coordinates
                const x = mol.getAtomX(i);
                const y = mol.getAtomY(i);
                const z = mol.getAtomZ(i) || 0;

                atoms.push({
                    id,
                    element: mol.getAtomLabel(i) || 'C', // Fallback
                    position: { x, y, z },
                    charge: mol.getAtomCharge(i)
                });
            }

            for (let i = 0; i < bondCount; i++) {
                const a1 = mol.getBondAtom(0, i);
                const a2 = mol.getBondAtom(1, i);
                const order = mol.getBondOrder(i);

                bonds.push({
                    id: nanoid(),
                    from: atomIdMap.get(a1)!,
                    to: atomIdMap.get(a2)!,
                    order: order > 3 ? 1 : order // Clamp aromatic/etc to 1 for now or handle simple types
                });
            }

            // Save to Supabase
            const { error } = await supabase
                .from('user_molecules')
                .insert({
                    user_id: userId,
                    name: file.name.split('.')[0],
                    atoms,
                    bonds,
                    smiles: mol.toSmiles(),
                    molfile: mol.toMolfile()
                });

            if (error) throw error;

            alert('Molecule imported successfully!');
            onImportSuccess();

        } catch (error) {
            console.error('Import failed:', error);
            alert('Failed to import molecule.');
        } finally {
            if (fileInputRef.current) fileInputRef.current.value = '';
        }
    };

    return (
        <>
            <input
                type="file"
                ref={fileInputRef}
                onChange={handleFileChange}
                accept=".mol,.sdf,.smi,.smiles,.txt"
                className="hidden"
            />
            <button
                onClick={() => {
                    if (!userId) {
                        alert('Please sign in to import molecules');
                        return;
                    }
                    fileInputRef.current?.click();
                }}
                className="btn-secondary flex items-center gap-2 px-3 py-1.5"
                title="Import Molecule (.mol, .sdf, .smiles)"
            >
                <Upload size={16} />
                <span>Import</span>
            </button>
        </>
    );
}
