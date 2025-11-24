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
 * Left sidebar tool buttons panel
 */
export default function ToolButtons() {
  const currentTool = useLabStore(s => s.currentTool)
  const setTool = useLabStore(s => s.setTool)

  return (
    <motion.div
      initial={{ opacity: 0, y: 10 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.4, ease: 'easeOut' }}
      className="flex flex-col items-center py-4 gap-3 w-full"
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
            onClick={(e) => {
              // Ripple effect
              const circle = document.createElement('span')
              circle.className = 'ripple'
              const rect = e.currentTarget.getBoundingClientRect()
              const size = Math.max(rect.width, rect.height)
              circle.style.width = circle.style.height = `${size}px`
              e.currentTarget.appendChild(circle)
              setTimeout(() => circle.remove(), 400)
              
              setTool(tool.id)
            }}
            className={`
              relative w-10 h-10 rounded-lg flex items-center justify-center border
              transition-colors overflow-hidden
              ${active
                ? 'border-[#4676ff] text-[#4676ff] bg-white shadow-sm'
                : 'border-[#ccc] text-gray-600'}
              hover:bg-white
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

