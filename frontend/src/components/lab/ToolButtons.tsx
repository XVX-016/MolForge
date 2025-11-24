import React from 'react'
import { motion } from 'framer-motion'
import { FiMousePointer, FiPlus, FiLink, FiTrash2 } from 'react-icons/fi'
import { useLabStore } from '../../store/labStore'
import type { ToolName } from '../../types/molecule'

const tools: Array<{ id: ToolName; icon: React.ReactNode; label: string }> = [
  { id: 'select', icon: <FiMousePointer size={18} />, label: 'Select' },
  { id: 'add_atom', icon: <FiPlus size={18} />, label: 'Add Atom' },
  { id: 'bond', icon: <FiLink size={18} />, label: 'Bond' },
  { id: 'delete', icon: <FiTrash2 size={18} />, label: 'Delete' },
]

/**
 * Floating vertical tool buttons panel (left side)
 * Soft embossed style matching target design
 */
export default function ToolButtons() {
  const currentTool = useLabStore(s => s.currentTool)
  const setTool = useLabStore(s => s.setTool)

  return (
    <motion.div
      initial={{ x: -20, opacity: 0 }}
      animate={{ x: 0, opacity: 1 }}
      transition={{ duration: 0.4, ease: 'easeOut' }}
      className="fixed top-1/2 left-10 -translate-y-1/2 z-20
        rounded-2xl p-4 bg-[#f5f1e8] border border-[#e3dfd6]
        shadow-[0_4px_20px_rgba(0,0,0,0.08)] flex flex-col gap-3"
    >
      {tools.map((tool, index) => {
        const active = currentTool === tool.id

        return (
          <motion.button
            key={tool.id}
            initial={{ opacity: 0, scale: 0.8 }}
            animate={{ opacity: 1, scale: 1 }}
            transition={{ duration: 0.3, delay: index * 0.05 }}
            whileHover={{ scale: 1.1 }}
            whileTap={{ scale: 0.95 }}
            onClick={() => setTool(tool.id)}
            className={`
              w-12 h-12 flex items-center justify-center rounded-full
              border transition-all duration-200
              ${active
                ? 'bg-white border-[#d4cec4] shadow-[inset_0_2px_4px_rgba(0,0,0,0.1)] text-[#4676ff]'
                : 'bg-[#faf8f3] border-[#e3dfd6] text-gray-600 hover:bg-white'}
            `}
            title={tool.label}
          >
            {tool.icon}
          </motion.button>
        )
      })}
    </motion.div>
  )
}

