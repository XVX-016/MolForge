
import { useEffect } from 'react';
import { motion } from 'framer-motion';
import { useStudioStore } from '../store/studioStore';
import { useHistoryStore } from '../store/historyStore';
import { BENZENE } from '../utils/defaultMolecules';
import ModeSwitcher from '../components/studio/ModeSwitcher';
import ActionsPanel from '../components/studio/ActionsPanel';
import MolecularWorkspace from '../components/studio/MolecularWorkspace';
import PropertyPanel from '../components/studio/PropertyPanel';

export default function StudioPage() {
  const {
    isCanvasInitialized,
    setCanvasInitialized
  } = useStudioStore();

  const { present, init, undo, redo } = useHistoryStore();

  // One-time initialization of default molecule (Benzene)
  useEffect(() => {
    if (!isCanvasInitialized && !present) {
      console.log('Initializing Studio with Default Molecule (Benzene)');
      init(BENZENE);
      setCanvasInitialized(true);
    }
  }, [isCanvasInitialized, present, init, setCanvasInitialized]);

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
    <motion.div
      initial={{ opacity: 0 }}
      animate={{ opacity: 1 }}
      transition={{ duration: 0.3 }}
      className="space-y-6 max-w-[1700px] mx-auto p-6 h-[calc(100vh-4rem)] flex flex-col"
    >
      <div className="flex flex-col md:flex-row items-center justify-between gap-4 shrink-0">
        <div>
          <h1 className="text-3xl font-extrabold text-black tracking-tight">MolForge Studio</h1>
          <p className="text-darkGrey text-sm">Design • Optimize • Simulate molecular systems</p>
        </div>

        <ModeSwitcher />
      </div>

      <div className="flex-1 grid grid-cols-1 lg:grid-cols-12 gap-8 min-h-0 px-4">

        {/* LEFT SIDEBAR: Actions & Properties */}
        <div className="lg:col-span-3 flex flex-col gap-6 h-full min-h-0">
          <div className="flex-1 min-h-0 rounded-2xl overflow-hidden shadow-sm border border-lightGrey">
            <ActionsPanel />
          </div>
          <div className="shrink-0 h-[280px]">
            <PropertyPanel />
          </div>
        </div>

        {/* MAIN WORKSPACE: More space, less sidebar dominance */}
        <div className="lg:col-span-9 h-full min-h-0 relative bg-gray-50/30 rounded-3xl border border-lightGrey/30 overflow-hidden">
          <MolecularWorkspace />
        </div>

      </div>
    </motion.div>
  );
}
