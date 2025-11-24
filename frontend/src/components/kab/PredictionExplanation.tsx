import React, { useState, useEffect } from 'react'
import { MoleculeGraph } from '@biosynth/engine'

interface PredictionExplanationProps {
  molecule: MoleculeGraph
  type: 'binding' | 'analysis' | 'knowledge'
  siteId?: number
  alertId?: number
  ruleId?: string
  onClose: () => void
}

export default function PredictionExplanation({
  molecule,
  type,
  siteId,
  alertId,
  ruleId,
  onClose
}: PredictionExplanationProps) {
  const [explanation, setExplanation] = useState<string>('')
  const [loading, setLoading] = useState(true)

  useEffect(() => {
    const fetchExplanation = async () => {
      try {
        const atoms = Array.from(molecule.atoms.values()).map((atom, idx) => ({
          element: atom.element,
          position: atom.position,
          id: idx
        }))
        
        const bonds: any[] = []
        molecule.bonds.forEach((bond) => {
          const a1Idx = atoms.findIndex(a => a.id === bond.a1)
          const a2Idx = atoms.findIndex(a => a.id === bond.a2)
          if (a1Idx >= 0 && a2Idx >= 0) {
            bonds.push({
              a1: a1Idx,
              a2: a2Idx,
              order: bond.order
            })
          }
        })

        const response = await fetch('/api/kab/explain', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            molecule: { atoms, bonds },
            site_id: siteId,
            alert_id: alertId,
            rule_id: ruleId
          })
        })

        if (!response.ok) {
          throw new Error('Failed to fetch explanation')
        }

        const data = await response.json()
        setExplanation(data.explanation || 'No explanation available.')
      } catch (err) {
        setExplanation(err instanceof Error ? err.message : 'Failed to generate explanation.')
      } finally {
        setLoading(false)
      }
    }

    fetchExplanation()
  }, [molecule, siteId, alertId, ruleId])

  return (
    <div className="fixed inset-0 bg-black bg-opacity-50 flex items-center justify-center z-50" onClick={onClose}>
      <div
        className="bg-white rounded-lg shadow-xl max-w-2xl w-full mx-4 max-h-[80vh] overflow-y-auto"
        onClick={(e) => e.stopPropagation()}
      >
        <div className="p-6">
          <div className="flex items-center justify-between mb-4">
            <h3 className="text-lg font-semibold text-gray-900">Explanation</h3>
            <button
              onClick={onClose}
              className="text-gray-400 hover:text-gray-600 transition-colors"
            >
              <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
              </svg>
            </button>
          </div>

          {loading ? (
            <div className="flex items-center justify-center py-8">
              <div className="w-6 h-6 border-2 border-blue-600 border-t-transparent rounded-full animate-spin"></div>
              <span className="ml-2 text-sm text-gray-600">Generating explanation...</span>
            </div>
          ) : (
            <div className="prose max-w-none">
              <p className="text-gray-700 whitespace-pre-wrap">{explanation}</p>
            </div>
          )}
        </div>
      </div>
    </div>
  )
}

