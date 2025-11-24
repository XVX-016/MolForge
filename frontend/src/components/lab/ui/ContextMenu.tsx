import React from 'react'
import { motion, AnimatePresence } from 'framer-motion'

interface ContextMenuProps {
  visible: boolean
  x: number
  y: number
  atomId: string | null
  bondId: string | null
  onDeleteAtom?: (atomId: string) => void
  onDeleteBond?: (bondId: string) => void
  onClose: () => void
}

/**
 * Context menu component for right-click actions
 */
export default function ContextMenu({
  visible,
  x,
  y,
  atomId,
  bondId,
  onDeleteAtom,
  onDeleteBond,
  onClose,
}: ContextMenuProps) {
  if (!visible) return null

  return (
    <AnimatePresence>
      {visible && (
        <motion.div
          initial={{ opacity: 0, scale: 0.9 }}
          animate={{ opacity: 1, scale: 1 }}
          exit={{ opacity: 0, scale: 0.9 }}
          className="fixed z-50 bg-white rounded-lg shadow-lg border border-gray-200 py-1 min-w-[160px]"
          style={{ left: x, top: y }}
          onClick={(e) => e.stopPropagation()}
        >
          {atomId && onDeleteAtom && (
            <button
              onClick={() => {
                onDeleteAtom(atomId)
                onClose()
              }}
              className="w-full px-4 py-2 text-left text-sm text-gray-700 hover:bg-gray-100 transition-colors"
            >
              Delete Atom
            </button>
          )}
          {bondId && onDeleteBond && (
            <button
              onClick={() => {
                onDeleteBond(bondId)
                onClose()
              }}
              className="w-full px-4 py-2 text-left text-sm text-gray-700 hover:bg-gray-100 transition-colors"
            >
              Delete Bond
            </button>
          )}
          {atomId && (
            <div className="border-t border-gray-200 my-1" />
          )}
          <button
            onClick={onClose}
            className="w-full px-4 py-2 text-left text-sm text-gray-500 hover:bg-gray-100 transition-colors"
          >
            Cancel
          </button>
        </motion.div>
      )}
    </AnimatePresence>
  )
}

