
import { useEffect } from 'react';
import { motion } from 'framer-motion';
import { useStudioStore } from '../store/studioStore';
import ChatInterface from '../components/studio/ChatInterface';
import Studio3DScene from '../components/studio/Studio3DScene';
import Card from '../components/ui/Card';
import { BENZENE } from '../utils/defaultMolecules';
import { MoleculeGraph } from '@biosynth/engine';

// Helper to hydrate Benzene into a Graph
const createBenzeneGraph = () => {
  const mol = new MoleculeGraph();
  const idMap = new Map<string, string>();

  // Add atoms
  BENZENE.atoms.forEach((atom) => {
    const newId = mol.addAtom({
      element: atom.element as any,
      position: [atom.x, atom.y, atom.z]
    });
    idMap.set(atom.id, newId);
  });

  // Add bonds
  BENZENE.bonds.forEach((bond) => {
    const fromId = idMap.get(bond.a);
    const toId = idMap.get(bond.b);
    if (fromId && toId) {
      mol.addBond(fromId, toId, bond.order);
    }
  });

  return mol;
};

export default function StudioPage() {
  const {
    mode,
    setMode,
    molecule,
    setMolecule,
    isCanvasInitialized,
    setCanvasInitialized
  } = useStudioStore();

  // One-time initialization of default molecule (Benzene)
  useEffect(() => {
    if (!isCanvasInitialized && !molecule) {
      console.log('Initializing Studio with Default Molecule (Benzene)');
      const benzene = createBenzeneGraph();
      setMolecule(benzene);
      setCanvasInitialized(true);
    }
  }, [isCanvasInitialized, molecule, setMolecule, setCanvasInitialized]);

  return (
    <motion.div
      initial={{ opacity: 0 }}
      animate={{ opacity: 1 }}
      transition={{ duration: 0.3 }}
      className="space-y-6 max-w-[1600px] mx-auto p-6 h-[calc(100vh-4rem)] flex flex-col"
    >
      {/* Studio Header & Mode Switcher */}
      <div className="flex flex-col md:flex-row items-center justify-between gap-4 shrink-0">
        <div>
          <h1 className="text-3xl font-extrabold text-black tracking-tight">MolForge Studio</h1>
          <p className="text-darkGrey text-sm">Design • Optimize • Simulate molecular systems</p>
        </div>

        {/* Mode Switcher */}
        <div className="bg-white p-1.5 rounded-xl border border-lightGrey shadow-sm flex items-center gap-1">
          {(['design', 'optimize', 'simulate'] as const).map((m) => (
            <button
              key={m}
              onClick={() => setMode(m)}
              className={`px-6 py-2 rounded-lg text-sm font-bold transition-all ${mode === m
                  ? 'bg-black text-white shadow-md transform scale-105'
                  : 'text-midGrey hover:text-black hover:bg-offwhite'
                }`}
            >
              {m.charAt(0).toUpperCase() + m.slice(1)}
            </button>
          ))}
        </div>
      </div>

      {/* Main Workspace Layout */}
      <div className="flex-1 grid grid-cols-1 lg:grid-cols-12 gap-6 min-h-0">

        {/* Command Panel (Left) */}
        <div className="lg:col-span-3 flex flex-col h-full min-h-0">
          <Card className="flex-1 flex flex-col bg-white border-lightGrey shadow-sm overflow-hidden">
            <div className="p-4 border-b border-lightGrey bg-offwhite/50 backdrop-blur-sm">
              <div className="flex items-center justify-between">
                <span className="text-xs font-bold uppercase tracking-wider text-midGrey">Command Center</span>
                <span className={`px-2 py-0.5 rounded text-[10px] font-bold uppercase ${mode === 'design' ? 'bg-blue-100 text-blue-700' :
                    mode === 'optimize' ? 'bg-purple-100 text-purple-700' :
                      'bg-orange-100 text-orange-700'
                  }`}>
                  {mode} Mode
                </span>
              </div>
            </div>

            <div className="flex-1 overflow-hidden">
              <ChatInterface />
            </div>
          </Card>
        </div>

        {/* Studio Canvas (Right) - The Hero */}
        <div className="lg:col-span-9 h-full min-h-0 relative group">
          <Card className="h-full w-full overflow-hidden bg-gradient-to-br from-offwhite to-white border-lightGrey relative p-0">
            {/* Canvas Area */}
            <div className="absolute inset-0 z-0">
              <Studio3DScene
                mode={mode}
                molecule={molecule}
              />
            </div>

            {/* Canvas Overlays based on Mode */}
            <div className="absolute top-4 left-4 z-10 pointer-events-none">
              <div className="inline-flex items-center gap-2 px-3 py-1.5 bg-white/90 backdrop-blur border border-lightGrey rounded-lg shadow-sm">
                <div className={`w-2 h-2 rounded-full ${mode === 'design' ? 'bg-green-500 animate-pulse' : 'bg-midGrey'
                  }`} />
                <span className="text-xs font-medium text-darkGrey">
                  {mode === 'design' ? 'Live Editing Active' :
                    mode === 'optimize' ? 'Optimization Locked' :
                      'Simulation View'}
                </span>
              </div>
            </div>

          </Card>
        </div>

      </div>
    </motion.div>
  );
}
