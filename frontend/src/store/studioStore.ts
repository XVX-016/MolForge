
import { create } from 'zustand';
import { MoleculeGraph } from '@biosynth/engine';
import type { StudioMode, StudioMessage } from '../types/studio';

interface StudioState {
    mode: StudioMode;
    molecule: MoleculeGraph | null;
    messages: StudioMessage[];
    isCanvasInitialized: boolean;

    // Actions
    setMode: (mode: StudioMode) => void;
    setMolecule: (molecule: MoleculeGraph | null) => void;
    addMessage: (message: StudioMessage) => void;
    resetMessages: () => void;
    setCanvasInitialized: (initialized: boolean) => void;
}

export const useStudioStore = create<StudioState>((set) => ({
    mode: 'design',
    molecule: null,
    messages: [],
    isCanvasInitialized: false,

    setMode: (mode) => set({ mode }),
    setMolecule: (molecule) => set({ molecule }),
    addMessage: (message) => set((state) => ({ messages: [...state.messages, message] })),
    resetMessages: () => set({ messages: [] }),
    setCanvasInitialized: (initialized) => set({ isCanvasInitialized: initialized }),
}));
