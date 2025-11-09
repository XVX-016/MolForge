import { create } from 'zustand'
import { MoleculeGraph, MoleculeSerializer } from '@biosynth/engine'
import { UndoStack } from '@biosynth/engine'
import { useMoleculeStore } from './moleculeStore'

interface HistoryState {
  undoStack: UndoStack
  canUndo: boolean
  canRedo: boolean
}

const createUndoStack = () => {
  const stack = new UndoStack()
  return stack
}

export const useHistoryStore = create<HistoryState>(() => ({
  undoStack: createUndoStack(),
  canUndo: false,
  canRedo: false,
}))

/**
 * Push current molecule state to history
 */
export function pushState(): void {
  const store = useMoleculeStore.getState()
  const history = useHistoryStore.getState()
  
  if (store.currentMolecule) {
    history.undoStack.push(store.currentMolecule)
    useHistoryStore.setState({
      canUndo: history.undoStack.canUndo(),
      canRedo: history.undoStack.canRedo(),
    })
  }
}

/**
 * Undo last action
 */
export function undo(): void {
  const history = useHistoryStore.getState()
  const restored = history.undoStack.undo()
  
  if (restored) {
    useMoleculeStore.getState().setMolecule(restored)
    useHistoryStore.setState({
      canUndo: history.undoStack.canUndo(),
      canRedo: history.undoStack.canRedo(),
    })
  }
}

/**
 * Redo last undone action
 */
export function redo(): void {
  const history = useHistoryStore.getState()
  const restored = history.undoStack.redo()
  
  if (restored) {
    useMoleculeStore.getState().setMolecule(restored)
    useHistoryStore.setState({
      canUndo: history.undoStack.canUndo(),
      canRedo: history.undoStack.canRedo(),
    })
  }
}

/**
 * Clear history
 */
export function clearHistory(): void {
  const history = useHistoryStore.getState()
  history.undoStack.clear()
  useHistoryStore.setState({
    canUndo: false,
    canRedo: false,
  })
}

