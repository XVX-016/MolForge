import React from 'react'
import { motion, AnimatePresence } from 'framer-motion'
import type { Element } from '@biosynth/engine'
import { useMoleculeStore } from '../store/moleculeStore'

// Professional metallic element palette - matches AtomMesh colors
const ELEMENTS: Array<{ element: Element; label: string; baseColor: string; accentColor: string }> = [
  { element: 'H', label: 'H', baseColor: '#F6F7F8', accentColor: '#FFFFFF' }, // Soft white/ivory
  { element: 'C', label: 'C', baseColor: '#9DA3AE', accentColor: '#2B2E33' }, // Charcoal grey
  { element: 'N', label: 'N', baseColor: '#8B95A1', accentColor: '#5A6B7A' }, // Colder grey
  { element: 'O', label: 'O', baseColor: '#6B8FA3', accentColor: '#4A6FA5' }, // Soft desaturated blue
  { element: 'F', label: 'F', baseColor: '#A8B5B8', accentColor: '#7A8B8F' }, // Soft beige-grey
  { element: 'S', label: 'S', baseColor: '#C4B5A0', accentColor: '#B5A085' }, // Soft beige tint
  { element: 'P', label: 'P', baseColor: '#A8A39D', accentColor: '#8B7A6B' }, // Warm grey
  { element: 'Cl', label: 'Cl', baseColor: '#9FA8B3', accentColor: '#6B7A8B' }, // Light steel grey
  { element: 'Br', label: 'Br', baseColor: '#8B8B8B', accentColor: '#6B6B6B' }, // Medium grey
  { element: 'I', label: 'I', baseColor: '#7A7A7A', accentColor: '#5A5A5A' }, // Darker grey
]

export default function AtomPalette() {
  const tool = useMoleculeStore((state) => state.tool)
  const atomToAdd = useMoleculeStore((state) => state.atomToAdd)
  const setAtomToAdd = useMoleculeStore((state) => state.setAtomToAdd)

  if (tool !== 'add-atom') return null

  return (
    <AnimatePresence>
      <motion.div
        initial={{ opacity: 0, y: -20 }}
        animate={{ opacity: 1, y: 0 }}
        exit={{ opacity: 0, y: -20 }}
        className="fixed left-20 top-4 z-50 frosted-glass rounded-lg shadow-glass border border-chrome/20 p-3 backdrop-blur-md"
      >
        <div className="text-xs text-gray-700 mb-3 font-semibold">Select Element</div>
        <div className="grid grid-cols-5 gap-2">
          {ELEMENTS.map((el) => (
            <motion.button
              key={el.element}
              whileHover={{ scale: 1.08, y: -2 }}
              whileTap={{ scale: 0.95 }}
              onClick={() => setAtomToAdd(el.element)}
              className={`w-12 h-12 rounded-lg font-bold text-sm transition-all border-2 ${
                atomToAdd === el.element
                  ? 'ring-2 ring-indigo-500 ring-offset-1 shadow-lg border-indigo-400'
                  : 'border-gray-300 hover:border-gray-400 shadow-sm'
              }`}
              style={{ 
                backgroundColor: el.baseColor,
                color: atomToAdd === el.element ? '#000000' : '#2B2E33',
                boxShadow: atomToAdd === el.element ? `0 0 0 2px ${el.accentColor}` : undefined
              }}
              title={el.label}
            >
              {el.label}
            </motion.button>
          ))}
        </div>
      </motion.div>
    </AnimatePresence>
  )
}

