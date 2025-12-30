
import { create } from 'zustand';
import type { StudioMode, StudioMessage } from '../types/studio';
import { processAICommand } from '../lib/studio/geminiControl';
import { validateAction } from '../lib/studio/validator';
import { useHistoryStore } from './historyStore';
import * as mutations from '../lib/mutations';
import { apiClient } from '../api/api';
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

    // Actions
    setMode: (mode: StudioMode) => void;
    addMessage: (message: StudioMessage) => void;
    resetMessages: () => void;
    setCanvasInitialized: (initialized: boolean) => void;
    setCommandRunning: (running: boolean) => void;
    setSelection: (type: SelectionType, id: string | null) => void;
    runCommand: (input: string) => Promise<void>;
}

export const useStudioStore = create<StudioState>((set, get) => ({
    mode: 'design',
    messages: [],
    isCanvasInitialized: false,
    isCommandRunning: false,
    selection: { type: null, id: null },

    setMode: (mode) => set({ mode }),
    addMessage: (message) => set((state) => ({ messages: [...state.messages, message] })),
    resetMessages: () => set({ messages: [] }),
    setCanvasInitialized: (initialized) => set({ isCanvasInitialized: initialized }),
    setCommandRunning: (running) => set({ isCommandRunning: running }),
    setSelection: (type, id) => set({ selection: { type, id } }),

    runCommand: async (input: string) => {
        const { mode, setCommandRunning } = get();
        const history = useHistoryStore.getState();
        const currentMolecule = history.present || { atoms: [], bonds: [] };

        setCommandRunning(true);

        try {
            // 1. Process via Gemini (mocked/impl in geminiControl)
            const action = await processAICommand(input, currentMolecule, mode);

            // 2. Validate against Contract & Mode
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
                case 'NO_OP':
                    description = `No Action: ${action.reason}`;
                    break;
            }

            history.applyMutation(nextGraph, description, 'ai', mode);

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
            setCommandRunning(false);
        }
    }
}));
