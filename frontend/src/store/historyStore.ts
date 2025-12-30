
import { create } from 'zustand';
import type { MoleculeGraph } from '../types/molecule';

export interface ActionLogEntry {
  id: string;
  timestamp: number;
  mode: "design" | "optimize" | "simulate";
  description: string;
  source: "user" | "ai" | "system";
}

interface HistoryState {
  past: MoleculeGraph[];
  present: MoleculeGraph | null;
  future: MoleculeGraph[];
  logs: ActionLogEntry[];

  // Core Actions
  init: (mol: MoleculeGraph) => void;
  applyMutation: (nextGraph: MoleculeGraph, description: string, source: ActionLogEntry['source'], mode?: ActionLogEntry['mode']) => void;
  undo: () => void;
  redo: () => void;
}

export const useHistoryStore = create<HistoryState>((set, get) => ({
  past: [],
  present: null,
  future: [],
  logs: [],

  init: (mol: MoleculeGraph) => {
    set({
      past: [],
      present: mol,
      future: [],
      logs: [{
        id: Date.now().toString(),
        timestamp: Date.now(),
        mode: 'design',
        description: 'Initialized Workspace',
        source: 'system'
      }]
    });
  },

  applyMutation: (nextGraph: MoleculeGraph, description: string, source: ActionLogEntry['source'], mode?: ActionLogEntry['mode']) => {
    const { present, past, logs } = get();

    // Create Log Entry
    const newLog: ActionLogEntry = {
      id: Date.now().toString(),
      timestamp: Date.now(),
      mode: mode || 'design',
      description: description,
      source: source
    };

    if (present) {
      set({
        past: [...past, present],
        present: nextGraph,
        future: [],
        logs: [...logs, newLog]
      });
    } else {
      set({
        present: nextGraph,
        past: [],
        future: [],
        logs: [...logs, newLog]
      });
    }
  },

  undo: () => {
    const { past, present, future, logs } = get();
    if (past.length === 0 || !present) return;

    const previous = past[past.length - 1];
    const newPast = past.slice(0, -1);

    const newLog: ActionLogEntry = {
      id: Date.now().toString(),
      timestamp: Date.now(),
      mode: 'design',
      description: 'Undid last action',
      source: 'user'
    };

    set({
      past: newPast,
      present: previous,
      future: [present, ...future],
      logs: [...logs, newLog]
    });
  },

  redo: () => {
    const { past, present, future, logs } = get();
    if (future.length === 0 || !present) return;

    const next = future[0];
    const newFuture = future.slice(1);

    const newLog: ActionLogEntry = {
      id: Date.now().toString(),
      timestamp: Date.now(),
      mode: 'design',
      description: 'Redid action',
      source: 'user'
    };

    set({
      past: [...past, present],
      present: next,
      future: newFuture,
      logs: [...logs, newLog]
    });
  }
}));
