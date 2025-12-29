
import { useEffect } from 'react';
import { motion } from 'framer-motion';
import { useStudioStore } from '../store/studioStore';
import { useHistoryStore } from '../store/historyStore';
import { BENZENE } from '../utils/defaultMolecules';
import ModeSwitcher from '../components/studio/ModeSwitcher';
import ActionsPanel from '../components/studio/ActionsPanel';
import MolecularWorkspace from '../components/studio/MolecularWorkspace';
import PropertyPanel from '../components/studio/PropertyPanel';
import { useStudioLayout } from '../lib/studio/layout';

export default function StudioPage() {
  const {
    mode,
    isCanvasInitialized,
    setCanvasInitialized
  } = useStudioStore();

  const { present, init, undo, redo } = useHistoryStore();
  const layout = useStudioLayout(mode);

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
      className="max-w-[1800px] mx-auto p-6 h-[calc(100vh-4rem)] flex flex-col gap-6"
    >
      {/* Header Area */}
      <div className="flex items-center justify-between shrink-0 px-4">
        <div>
          <h1 className="text-3xl font-black text-black tracking-tight uppercase">MolForge <span className="text-gray-400">Studio</span></h1>
        </div>
        <ModeSwitcher />
      </div>

      {/* Main Structural Layout */}
      <div className="flex-1 grid grid-cols-12 gap-8 min-h-0 px-4">

        {/* SIDEBAR: Conditionally rendered based on layout contract */}
        {layout.showActionsPanel && (
          <div className="col-span-3 flex flex-col gap-6 h-full min-h-0">
            <div className="flex-1 min-h-0 rounded-2xl overflow-hidden shadow-sm border border-lightGrey bg-white">
              <ActionsPanel />
            </div>

            {layout.showPropertyPanel && (
              <div className="shrink-0 h-[280px]">
                <PropertyPanel />
              </div>
            )}
          </div>
        )}

        {/* WORKSPACE: Takes remaining space */}
        <div className={`${layout.showActionsPanel ? 'col-span-9' : 'col-span-12'} h-full min-h-0 relative rounded-3xl overflow-hidden border border-lightGrey/30 shadow-inner bg-gray-50/20`}>
          <MolecularWorkspace />
        </div>

      </div>
    </motion.div>
  );
}
