/**
 * selectionManager.ts
 * Simple selection manager for the kernel. Not DOM-dependent.
 * Emits events: 'select', 'deselect', 'change'
 * 
 * This is a kernel-level selection manager, separate from the UI SelectionManager
 * in components/r3f/SelectionManager.ts which handles hover/drag state.
 */

type SelectionChangePayload = { selected: string | null };
type SelectionEventMap = {
  select: string;
  deselect: string;
  change: SelectionChangePayload;
};
type SelectionEventName = keyof SelectionEventMap;
type Handler<T> = (payload: T) => void;

export class KernelSelectionManager {
  private selectedAtomId: string | null = null;
  private listeners: {
    [K in SelectionEventName]: Set<Handler<SelectionEventMap[K]>>;
  } = {
    select: new Set(),
    deselect: new Set(),
    change: new Set(),
  };

  constructor() {}

  getSelectedAtomId(): string | null {
    return this.selectedAtomId;
  }

  selectAtom(id: string): void {
    if (this.selectedAtomId === id) return;
    this.selectedAtomId = id;
    this.emit('select', id);
    this.emit('change', { selected: id });
  }

  deselect(): void {
    if (!this.selectedAtomId) return;
    const prev = this.selectedAtomId;
    this.selectedAtomId = null;
    this.emit('deselect', prev);
    this.emit('change', { selected: null });
  }

  on<K extends SelectionEventName>(event: K, cb: Handler<SelectionEventMap[K]>): () => void {
    const set = this.listeners[event];
    set?.add(cb);
    return () => set?.delete(cb);
  }

  private emit<K extends SelectionEventName>(event: K, payload: SelectionEventMap[K]): void {
    const set = this.listeners[event];
    if (!set) return;
    for (const h of Array.from(set)) {
      try {
        h(payload);
      } catch (error) {
        // swallow handler errors (kernels shouldn't crash tests)
        console.warn('KernelSelectionManager handler error', error);
      }
    }
  }

  reset(): void {
    this.selectedAtomId = null;
    this.listeners.select.clear();
    this.listeners.deselect.clear();
    this.listeners.change.clear();
  }
}

// Export singleton instance
export const kernelSelectionManager = new KernelSelectionManager();

// Also export class for testing
export default KernelSelectionManager;

