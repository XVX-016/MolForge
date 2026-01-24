
import { create } from 'zustand';
import type { StudioMode, StudioMessage } from '../types/studio';
import { processAICommand } from '../lib/studio/geminiControl';
import { validateAction } from '../lib/studio/validator';
import { useHistoryStore } from './historyStore';
import * as mutations from '../lib/mutations';
import { apiClient } from '../api/api';
import { toSMILES, getMoleculeProperties, commitVersion, analyzeMolecule } from '../lib/api';
import type { MoleculeProperties, AnalysisResponse } from '../lib/api';
import type { MoleculeGraph } from '../types/molecule';

export type SelectionType = 'atom' | 'bond' | null;

interface StudioSelection {
    type: SelectionType;
    id: string | null;
}

interface StudioState {
    mode: StudioMode;
    messages: StudioMessage[];
    isCanvasInitialized: boolean;
    isCommandRunning: boolean;
    selection: StudioSelection;
    properties: MoleculeProperties | null;
    analysis: AnalysisResponse | null;
    isComputingProperties: boolean;
    activeVersionId: string | null;
    activeVersionIndex: number | null;
    isSyncing: boolean;
    isDirty: boolean;
    canCommit: boolean;

    // Actions
    setMode: (mode: StudioMode) => void;
    addMessage: (message: StudioMessage) => void;
    resetMessages: () => void;
    setCanvasInitialized: (initialized: boolean) => void;
    setCommandRunning: (running: boolean) => void;
    setSelection: (type: SelectionType, id: string | null) => void;
    runCommand: (input: string) => Promise<void>;
    fetchMoleculeProperties: (graph: MoleculeGraph) => Promise<void>;
    fetchAnalysis: (smiles: string) => Promise<void>;
    commitDraft: () => Promise<void>;
}

export const useStudioStore = create<StudioState>((set, get) => ({
    mode: 'design',
    messages: [],
    isCanvasInitialized: false,
    isCommandRunning: false,
    selection: { type: null, id: null },
    properties: null,
    analysis: null,
    isComputingProperties: false,
    activeVersionId: null,
    activeVersionIndex: null,
    isSyncing: false,
    isDirty: false,
    canCommit: false,

    setMode: (mode) => set({ mode }),
    addMessage: (message) => set((state) => ({ messages: [...state.messages, message] })),
    resetMessages: () => set({ messages: [] }),
    setCanvasInitialized: (initialized) => set({ isCanvasInitialized: initialized }),
    setCommandRunning: (running) => set({ isCommandRunning: running }),
    setSelection: (type, id) => set({ selection: { type, id } }),

    fetchMoleculeProperties: async (graph: MoleculeGraph) => {
        set({ isComputingProperties: true });
        try {
            const { smiles } = await toSMILES(graph);
            if (smiles) {
                const props = await getMoleculeProperties(smiles);
                set({
                    properties: props,
                    canCommit: !!props.canonical_smiles && !!props.inchikey
                });
                await get().fetchAnalysis(smiles);
            }
        } catch (err) {
            console.error('Failed to fetch molecular properties:', err);
            // Fallback UX rule: Keep last known if reachable, or show null
            // For now, we just log.
        } finally {
            set({ isComputingProperties: false });
        }
    },

    fetchAnalysis: async (smiles: string) => {
        try {
            const analysis = await analyzeMolecule(smiles);
            set({ analysis });
        } catch (err) {
            console.error('Failed to fetch analysis:', err);
        }
    },

    runCommand: async (input: string) => {
        const { mode, setCommandRunning } = get();
        const history = useHistoryStore.getState();
        const currentMolecule = history.present || { atoms: [], bonds: [] };

        setCommandRunning(true);

        try {
            // 1. Process via Gemini (Command Manual Mode)
            const action = await processAICommand(input, currentMolecule, mode, get().analysis);

            // 2. Validate against Contract & Baseline Heuristics
            if (action.type === 'SELECT_OPTIMIZATION_RULE') {
                const ruleId = (action.payload as any)?.rule_id;
                const availableRule = get().analysis?.suggestions.find(s => s.id === ruleId);

                if (!availableRule) {
                    throw new Error(`AI attempted to hallucinate rule: ${ruleId}`);
                }

                // Propose the mutation found in the analysis
                // In Phase 4, "Applying" an AI rule creates an OPTIMIZED draft state.
                // For this demo, we can just trigger the existing availableRule logic
                // But we must follow the "ORCHESTRATOR" rule - AI doesn't know how to move atoms.
                // It just says "DO THE THING RDKIT SAID".
            }

            const validation = validateAction(action, currentMolecule, mode);

            if (!validation.valid) {
                history.applyMutation(
                    currentMolecule,
                    `Action Rejected: ${validation.error}`,
                    'system',
                    mode
                );
                return;
            }

            // 3. Apply Mutation
            let nextGraph: MoleculeGraph = currentMolecule;
            let description = action.reason || `Executed ${action.type}`;

            switch (action.type) {
                case 'CREATE_MOLECULE':
                    // If backend return full atoms/bonds, use them. 
                    // Backend might return simplified IDs, we might need to map them if user had something already.
                    // But in CREATE, we usually start fresh or replace entirely.
                    nextGraph = action.payload;
                    break;
                case 'ADD_ATOM':
                    nextGraph = mutations.addAtom(currentMolecule, action.payload.element, action.payload.position);
                    break;
                case 'REMOVE_ATOM':
                    nextGraph = mutations.deleteAtom(currentMolecule, action.payload.atomId);
                    break;
                case 'REPLACE_ATOM':
                    nextGraph = mutations.updateAtomElement(currentMolecule, action.payload.atomId, action.payload.newElement);
                    break;
                case 'ADD_BOND':
                    nextGraph = mutations.addBond(currentMolecule, action.payload.from, action.payload.to, action.payload.order as any);
                    break;
                case 'REMOVE_BOND':
                    nextGraph = mutations.deleteBond(currentMolecule, action.payload.bondId);
                    break;
                case 'OPTIMIZE_GEOMETRY':
                    try {
                        const response = await apiClient.post('/api/molecule/generate-3d', {
                            molecule: currentMolecule,
                            optimize: true
                        });
                        if (response.data.molecule) {
                            nextGraph = response.data.molecule;
                            description = "Geometry optimized via MMFF/UFF.";
                        }
                    } catch (optError: any) {
                        console.error('Optimization failed:', optError);
                        description = `Optimization failed: ${optError.message}`;
                    }
                    break;
                case 'SIMULATE_REACTION':
                    description = "Initiating reaction simulation...";
                    // Simulation logic would go here
                    break;
                case 'SELECT_OPTIMIZATION_RULE':
                    // This action types the AI design intent to a deterministic kernel suggest
                    // We can reuse history logic but mark it as 'ai' intent.
                    description = `AI Proposal: ${action.reason}`;
                    // In a full implementation, this might lead to a ghost/proposal state.
                    // For now, we apply it to demonstrate the "locked" orchestration.
                    break;
                case 'NO_OP':
                    description = `AI Response: ${action.reason}`;
                    break;
            }

            history.applyMutation(nextGraph, description, 'ai', mode);
            set({ isDirty: true });
            // Trigger property re-computation
            get().fetchMoleculeProperties(nextGraph);

        } catch (error: any) {
            console.error('runCommand Error:', error);

            let sanitizedReason = "The AI Architect encountered an internal processing error.";
            const errorMsg = error.message || "";

            if (errorMsg.includes('API key') || errorMsg.includes('400')) {
                sanitizedReason = "AI Orchestration services are currently restricted. Please check system configuration (API Key).";
            } else if (errorMsg.includes('timeout') || errorMsg.includes('504')) {
                sanitizedReason = "Connection to the molecular reasoning engine timed out. Re-submit your request.";
            } else if (errorMsg.includes('429')) {
                sanitizedReason = "Molecular reasoning throughput limit reached. Rate limiting active.";
            }

            history.applyMutation(
                currentMolecule,
                `Architect Unavailable: ${sanitizedReason}`,
                'system',
                mode
            );
        } finally {
            set({ isCommandRunning: false });
        }
    },

    commitDraft: async () => {
        const { activeVersionId } = get();
        const history = useHistoryStore.getState();
        const currentMolecule = history.present;

        if (!currentMolecule || currentMolecule.atoms.length === 0) return;

        set({ isSyncing: true });
        try {
            const { smiles } = await toSMILES(currentMolecule);
            const response = await commitVersion(smiles, currentMolecule, activeVersionId);

            set({
                activeVersionId: response.version_id,
                activeVersionIndex: response.version_index,
                properties: response.properties,
                isDirty: false
            });

            useHistoryStore.getState().applyMutation(
                currentMolecule,
                `Snapshot Committed: Version ${response.version_index}`,
                'system',
                get().mode
            );
        } catch (err) {
            console.error('Failed to commit molecular version:', err);
        } finally {
            set({ isSyncing: false });
        }
    }
}));
