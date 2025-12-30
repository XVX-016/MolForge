import { useEffect } from 'react';
import AIPanel from '../components/studio/AIPanel';
import MolecularWorkspace from '../components/studio/MolecularWorkspace';
import { useHistoryStore } from '../store/historyStore';

export default function StudioPage() {
  const { undo, redo } = useHistoryStore();

  // Keyboard Shortcuts for Undo/Redo
  useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      // Undo: Ctrl+Z or Cmd+Z
      if ((e.metaKey || e.ctrlKey) && e.key === 'z') {
        e.preventDefault();
        if (e.shiftKey) {
          redo();
        } else {
          undo();
        }
      }
      // Redo: Ctrl+Y or Cmd+Y
      if ((e.metaKey || e.ctrlKey) && e.key === 'y') {
        e.preventDefault();
        redo();
      }
    };

    window.addEventListener('keydown', handleKeyDown);
    return () => window.removeEventListener('keydown', handleKeyDown);
  }, [undo, redo]);

  return (
    <div className="flex h-full overflow-hidden">
      {/* AI Control Plane (Docked Left) */}
      <aside className="w-[380px] shrink-0 border-r border-[#1F2937] bg-[#0F1013]">
        <AIPanel />
      </aside>

      {/* The Stage (Canvas) */}
      <section className="flex-1 relative bg-black">
        <MolecularWorkspace />
      </section>
    </div>
  );
}
