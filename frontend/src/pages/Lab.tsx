import React, { useState, useEffect } from 'react';
import MoleculeViewer from '../components/MoleculeViewer';
import ToolPanel from '../components/ToolPanel';
import PropertiesPanel from '../components/PropertiesPanel';
import { TemplatePanel } from '../components/TemplatePanel';
import Button from '../components/ui/Button';
import { useMoleculeStore } from '../store/moleculeStore';
import { saveMolecule as saveMoleculeToSupabase } from '../lib/supabaseMoleculeStore';
import { moleculeToJSON, getCanvasThumbnail } from '../lib/engineAdapter';
import { MoleculeSerializer, MoleculeGraph } from '@biosynth/engine';
import { supabase } from '../supabase';
import { convertSMILESToMolfile } from '../lib/api';

export default function Lab() {
	const molecule = useMoleculeStore((s) => s.currentMolecule);
	const backendPredictions = useMoleculeStore((s) => s.backendPredictions);
	const [saving, setSaving] = useState(false);
	const [userId, setUserId] = useState<string | null>(null);

	useEffect(() => {
		if (!supabase) return;

		// Check current session
		supabase.auth.getSession().then(({ data: { session } }) => {
			if (session?.user) {
				setUserId(session.user.id);
			}
		});

		// Listen for auth changes
		const {
			data: { subscription },
		} = supabase.auth.onAuthStateChange((_event, session) => {
			if (session?.user) {
				setUserId(session.user.id);
			} else {
				setUserId(null);
			}
		});

		return () => subscription.unsubscribe();
	}, []);

	const handleSave = async () => {
		if (!molecule) return;
		
		if (!userId) {
			alert('Please sign in to save molecules. Visit /supabase-test to sign in.');
			return;
		}

		const name = prompt('Enter molecule name:');
		if (!name) return;
		
		setSaving(true);
		try {
			const jsonGraph = moleculeToJSON(molecule);
			const smiles = MoleculeSerializer.toSMILES(molecule);
			const thumbnail = getCanvasThumbnail();
			const properties = backendPredictions ? JSON.stringify(backendPredictions) : undefined;
			
			// Calculate formula if possible
			const formula = molecule ? calculateFormula(molecule) : undefined;
			
			// Generate molfile from SMILES if available
			let molfile: string | undefined = undefined;
			if (smiles && smiles.trim()) {
				try {
					const result = await convertSMILESToMolfile(smiles);
					molfile = result.molfile;
				} catch (e: any) {
					console.warn('Failed to generate molfile:', e);
					// Continue without molfile - it's optional
				}
			}
			
			await saveMoleculeToSupabase(userId, {
				name,
				smiles: smiles || undefined,
				formula,
				molfile,
				json_graph: jsonGraph,
				properties,
				thumbnail_b64: thumbnail || undefined,
			});
			alert(`Molecule "${name}" saved successfully!`);
		} catch (e: any) {
			console.error('Save error:', e);
			alert(`Failed to save molecule: ${e.message || 'Unknown error'}`);
		} finally {
			setSaving(false);
		}
	};

	// Helper function to calculate formula
	const calculateFormula = (mol: MoleculeGraph | null): string | undefined => {
		if (!mol) return undefined;
		const elementCounts: Record<string, number> = {};
		mol.atoms.forEach((atom) => {
			elementCounts[atom.element] = (elementCounts[atom.element] || 0) + 1;
		});
		return Object.entries(elementCounts)
			.map(([el, count]) => count > 1 ? `${el}${count}` : el)
			.join('');
	};

	return (
		<div className="lab-layout">
			<TemplatePanel />
			<div className="grid grid-cols-12 gap-4 flex-1">
				<div className="col-span-12 lg:col-span-3">
					<div className="bg-white rounded-xl shadow-neon border border-lightGrey p-4">
						<ToolPanel />
					</div>
				</div>
				<div className="col-span-12 lg:col-span-6">
					<div className="relative rounded-xl shadow-neon border border-lightGrey p-2 h-[70vh] lg:h-[78vh] bg-offwhite">
						<MoleculeViewer />
						<div className="absolute bottom-3 right-3">
							<Button onClick={handleSave} disabled={!molecule || saving || !userId}>
								{saving ? 'Saving...' : userId ? 'Save to Library' : 'Sign in to Save'}
							</Button>
						</div>
					</div>
				</div>
				<div className="col-span-12 lg:col-span-3">
					<div className="bg-white rounded-xl shadow-neon border border-lightGrey p-4">
						<PropertiesPanel />
					</div>
				</div>
			</div>
		</div>
	);
}


