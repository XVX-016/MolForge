
import { create } from 'zustand';
import type { StudioMode, StudioMessage } from '../types/studio';

export type SelectionType = 'atom' | 'bond' | null;

interface StudioSelection {
    type: SelectionType;
    id: string | null;
}

interface StudioState {
    mode: StudioMode;
    messages: StudioMessage[];
    isCanvasInitialized: boolean;
    selection: StudioSelection;

    // Actions
    setMode: (mode: StudioMode) => void;
    addMessage: (message: StudioMessage) => void;
    resetMessages: () => void;
    setCanvasInitialized: (initialized: boolean) => void;
    setSelection: (type: SelectionType, id: string | null) => void;
}

export const useStudioStore = create<StudioState>((set) => ({
    mode: 'design',
    messages: [],
    isCanvasInitialized: false,
    selection: { type: null, id: null },

    setMode: (mode) => set({ mode }),
    addMessage: (message) => set((state) => ({ messages: [...state.messages, message] })),
    resetMessages: () => set({ messages: [] }),
    setCanvasInitialized: (initialized) => set({ isCanvasInitialized: initialized }),
    setSelection: (type, id) => set({ selection: { type, id } }),
}));
