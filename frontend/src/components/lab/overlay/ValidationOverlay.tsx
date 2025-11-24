import React from 'react'
import { motion, AnimatePresence } from 'framer-motion'
import { useLabStore } from '../../../store/labStore'
import { validateStructure } from '../../../utils/chemistry/validateStructure'

/**
 * Validation overlay showing structure issues
 */
export default function ValidationOverlay() {
  const mol = useLabStore(s => s.molecule)
  const showValidation = useLabStore(s => s.showValidation ?? true)

  if (!showValidation) return null

  const issues = validateStructure(mol)

  if (issues.length === 0) return null

  return (
    <motion.div
      initial={{ opacity: 0, y: -10 }}
      animate={{ opacity: 1, y: 0 }}
      className="absolute top-20 right-4 bg-white rounded-lg shadow-lg border border-yellow-300 p-4 max-w-sm z-50"
    >
      <div className="flex items-center gap-2 mb-2">
        <span className="text-yellow-600">⚠</span>
        <h4 className="text-sm font-semibold text-gray-800">Structure Issues</h4>
      </div>
      <div className="space-y-1 text-xs">
        {issues.map((issue, i) => (
          <div
            key={i}
            className={`p-2 rounded ${
              issue.type === 'overbonded' ? 'bg-red-50 text-red-700' : 'bg-yellow-50 text-yellow-700'
            }`}
          >
            {issue.message}
          </div>
        ))}
      </div>
    </motion.div>
  )
}

