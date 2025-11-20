import React, { useState } from 'react'
import { motion } from 'framer-motion'
import { useInView } from 'react-intersection-observer'
import type { MoleculeItem } from '../lib/api'
import Molecule3DViewer from './Molecule3DViewer'

interface MoleculeCardProps {
  item: MoleculeItem & { molfile?: string | null; formula?: string | null }
  onOpen: () => void
  onDelete: () => void
}

export default function MoleculeCard({ item, onOpen, onDelete }: MoleculeCardProps) {
  const [loaded, setLoaded] = useState(false)
  
  // Lazy load 3D viewer only when card enters viewport
  const { ref, inView } = useInView({
    triggerOnce: true,
    threshold: 0.1,
    rootMargin: '50px', // Start loading slightly before card is visible
  })

  const formatDate = (dateStr: string) => {
    const date = new Date(dateStr)
    return date.toLocaleDateString('en-US', { month: 'short', day: 'numeric', year: 'numeric' })
  }

  const has3DData = item.molfile || item.smiles
  const showThumbnail = item.thumbnail_b64 && (!loaded || !has3DData)

  return (
    <motion.div
      ref={ref}
      layout
      initial={{ opacity: 0, y: 20 }}
      animate={{ opacity: 1, y: 0 }}
      whileHover={{ y: -6, boxShadow: '0px 10px 30px rgba(0,0,0,0.10)' }}
      className="bg-white p-4 rounded-xl shadow-neon border border-lightGrey hover:shadow-neon-hover hover:border-midGrey transition-all"
    >
      {/* 3D Viewer or Thumbnail with lazy loading */}
      <div className="h-40 bg-offwhite rounded-lg overflow-hidden mb-3 relative">
        {/* Thumbnail fallback (blurred background while loading) */}
        {showThumbnail && (
          <img
            src={item.thumbnail_b64}
            alt={item.name}
            className={`w-full h-full object-contain transition-opacity duration-300 ${
              loaded && has3DData ? 'opacity-30 blur-sm' : 'opacity-100'
            }`}
          />
        )}
        
        {/* Lazy-mounted 3D Viewer - only renders when in viewport */}
        {inView && has3DData && (
          <div className={loaded ? 'opacity-100' : 'opacity-0'}>
            <Molecule3DViewer
              molfile={item.molfile || undefined}
              smiles={item.smiles || undefined}
              height={160}
              backgroundColor="#f5f5f5"
              autorotate={false}
              onLoaded={() => setLoaded(true)}
            />
          </div>
        )}
        
        {/* Fallback when no data available */}
        {!has3DData && !showThumbnail && (
          <div className="w-full h-full flex items-center justify-center text-midGrey text-sm">
            No preview
          </div>
        )}
      </div>

      {/* Info */}
      <h3 className="font-semibold text-lg text-black mb-1 truncate">{item.name}</h3>
      {item.formula && (
        <div className="text-sm text-darkGrey mb-1">
          Formula: <span className="font-mono">{item.formula}</span>
        </div>
      )}
      {item.smiles && (
        <div className="text-xs text-midGrey mb-2 font-mono truncate" title={item.smiles}>
          SMILES: {item.smiles}
        </div>
      )}
      <div className="text-xs text-midGrey mb-3">{formatDate(item.created_at)}</div>

      {/* Actions */}
      <div className="flex gap-2">
        <motion.button
          whileHover={{ scale: 1.05 }}
          whileTap={{ scale: 0.95 }}
          onClick={onOpen}
          className="flex-1 px-3 py-2 rounded-lg bg-black text-white text-sm font-medium hover:bg-darkGrey hover:shadow-neon transition-all"
        >
          Open in Lab
        </motion.button>
        <motion.button
          whileHover={{ scale: 1.05 }}
          whileTap={{ scale: 0.95 }}
          onClick={onDelete}
          className="px-3 py-2 rounded-lg bg-white text-darkGrey hover:text-black hover:bg-offwhite border border-lightGrey text-sm font-medium transition-all"
        >
          Delete
        </motion.button>
      </div>
    </motion.div>
  )
}

