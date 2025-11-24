import React, { useState } from 'react'
import { motion } from 'framer-motion'
import Navbar from '../../components/Navbar'
import LabPanel from '../../components/lab/LabPanel'
import LabViewer from '../../components/lab/LabViewer'
import ToolButtons from '../../components/lab/ToolButtons'
import { useLabStore } from '../../store/labStore'

/**
 * Modern centered Lab layout matching target design:
 * - Floating tool buttons (left)
 * - Centered soft beige panel
 * - White 3D workspace
 * - Minimal controls
 */
export default function NewLabLayout() {
  const [mobileOpen, setMobileOpen] = useState(false)
  const currentTool = useLabStore(s => s.currentTool)
  const currentElement = useLabStore(s => s.currentElement)

  return (
    <div className="min-h-screen bg-white overflow-hidden">
      {/* Standard Navbar */}
      <Navbar onToggleMenu={() => setMobileOpen((v) => !v)} />

      {/* Floating Tool Buttons (left side) */}
      <ToolButtons />

      {/* Centered Workspace Panel */}
      <LabPanel>
        {/* Top Controls */}
        <motion.div
          initial={{ y: -10, opacity: 0 }}
          animate={{ y: 0, opacity: 1 }}
          transition={{ duration: 0.3 }}
          className="flex justify-between items-center mb-4"
        >
          <button
            className="px-4 py-2 rounded-xl bg-[#eae5dd] border border-[#d4cec4]
              text-sm font-medium text-gray-700 hover:bg-[#dfd9d0] transition-colors
              shadow-[inset_0_1px_2px_rgba(0,0,0,0.05)]"
            style={{ textTransform: 'uppercase' }}
          >
            {currentTool === 'add_atom' ? `Add ${currentElement}` : 'Add Atom'}
          </button>

          <div className="w-10 h-10 rounded-full border border-[#d4cec4] 
            bg-white flex items-center justify-center text-gray-500
            shadow-[inset_0_1px_2px_rgba(0,0,0,0.05)]">
            →
          </div>
        </motion.div>

        {/* Main 3D Viewer */}
        <motion.div
          initial={{ opacity: 0, scale: 0.98 }}
          animate={{ opacity: 1, scale: 1 }}
          transition={{ duration: 0.5 }}
          className="rounded-2xl overflow-hidden border border-[#d4cec4] bg-white"
          style={{ height: '500px' }}
        >
          <LabViewer />
        </motion.div>

        {/* Bottom Warning/Status */}
        <motion.div
          initial={{ y: 10, opacity: 0 }}
          animate={{ y: 0, opacity: 1 }}
          transition={{ duration: 0.3, delay: 0.2 }}
          className="mt-4 p-3 rounded-xl bg-[#eae5dd] border border-[#d4cec4] 
            flex items-center gap-2 text-sm text-gray-700
            shadow-[inset_0_1px_2px_rgba(0,0,0,0.05)]"
        >
          <span>⚠</span>
          <span>UNSTABLE STRUCTURE</span>
        </motion.div>
      </LabPanel>
    </div>
  )
}
