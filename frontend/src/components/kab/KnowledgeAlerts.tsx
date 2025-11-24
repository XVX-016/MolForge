import React, { useState, useEffect } from 'react'
import { MoleculeGraph } from '@biosynth/engine'
import { useMoleculeStore } from '../../store/moleculeStore'
import PredictionExplanation from './PredictionExplanation'

interface KnowledgeRule {
  rule: string
  description: string
  affected_atoms: number[]
  confidence: number
  category: string
}

interface MLPredictions {
  drug_likeness: number
  bioavailability: number
  toxicity_risk: number
  metabolic_stability: number
}

interface KnowledgeData {
  rules: KnowledgeRule[]
  ml_predictions: MLPredictions
  rule_count: number
}

interface KnowledgeAlertsProps {
  molecule: MoleculeGraph
  onHighlightAtoms?: (atomIndices: number[]) => void
}

export default function KnowledgeAlerts({ molecule, onHighlightAtoms }: KnowledgeAlertsProps) {
  const [data, setData] = useState<KnowledgeData | null>(null)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const [selectedRule, setSelectedRule] = useState<string | null>(null)
  const [showExplanation, setShowExplanation] = useState(false)
  const selectAtom = useMoleculeStore((s) => s.selectAtom)

  useEffect(() => {
    if (!molecule || molecule.atoms.size === 0) {
      setData(null)
      return
    }

    const fetchKnowledge = async () => {
      setLoading(true)
      setError(null)
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

        const response = await fetch('/api/kab/knowledge', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ atoms, bonds })
        })

        if (!response.ok) {
          throw new Error('Failed to fetch knowledge rules')
        }

        const result = await response.json()
        setData(result)
      } catch (err) {
        setError(err instanceof Error ? err.message : 'Unknown error')
        console.error('Knowledge fetch error:', err)
      } finally {
        setLoading(false)
      }
    }

    fetchKnowledge()
  }, [molecule])

  const getCategoryColor = (category: string) => {
    const colors: Record<string, string> = {
      reactivity: 'bg-blue-100 text-blue-800',
      stability: 'bg-green-100 text-green-800',
      drug_likeness: 'bg-purple-100 text-purple-800',
      stereochemistry: 'bg-yellow-100 text-yellow-800'
    }
    return colors[category] || 'bg-gray-100 text-gray-800'
  }

  const getConfidenceColor = (confidence: number) => {
    if (confidence >= 0.8) return 'text-green-600'
    if (confidence >= 0.6) return 'text-yellow-600'
    return 'text-red-600'
  }

  if (loading) {
    return (
      <div className="flex items-center justify-center py-8">
        <div className="w-6 h-6 border-2 border-blue-600 border-t-transparent rounded-full animate-spin"></div>
        <span className="ml-2 text-sm text-gray-600">Applying knowledge rules...</span>
      </div>
    )
  }

  if (error) {
    return (
      <div className="p-4 bg-red-50 border border-red-200 rounded-lg">
        <p className="text-sm text-red-800">Error: {error}</p>
      </div>
    )
  }

  if (!data) {
    return (
      <div className="p-4 text-sm text-gray-500 text-center">
        No knowledge data available.
      </div>
    )
  }

  return (
    <div className="space-y-6">
      {/* ML Predictions */}
      {data.ml_predictions && (
        <div>
          <h4 className="text-sm font-semibold text-gray-900 mb-3">ML Predictions</h4>
          <div className="grid grid-cols-2 gap-3">
            <div className="p-3 bg-purple-50 rounded-lg">
              <div className="text-xs text-gray-600 mb-1">Drug Likeness</div>
              <div className="text-lg font-semibold">{(data.ml_predictions.drug_likeness * 100).toFixed(1)}%</div>
            </div>
            <div className="p-3 bg-blue-50 rounded-lg">
              <div className="text-xs text-gray-600 mb-1">Bioavailability</div>
              <div className="text-lg font-semibold">{(data.ml_predictions.bioavailability * 100).toFixed(1)}%</div>
            </div>
            <div className="p-3 bg-red-50 rounded-lg">
              <div className="text-xs text-gray-600 mb-1">Toxicity Risk</div>
              <div className="text-lg font-semibold">{(data.ml_predictions.toxicity_risk * 100).toFixed(1)}%</div>
            </div>
            <div className="p-3 bg-green-50 rounded-lg">
              <div className="text-xs text-gray-600 mb-1">Metabolic Stability</div>
              <div className="text-lg font-semibold">{(data.ml_predictions.metabolic_stability * 100).toFixed(1)}%</div>
            </div>
          </div>
        </div>
      )}

      {/* Knowledge Rules */}
      {data.rules && data.rules.length > 0 && (
        <div>
          <h4 className="text-sm font-semibold text-gray-900 mb-3">Knowledge Rules ({data.rules.length})</h4>
          <div className="space-y-2">
            {data.rules.map((rule, idx) => (
              <div
                key={idx}
                className="p-3 border border-gray-200 rounded-lg hover:border-gray-300 hover:bg-gray-50 transition-all cursor-pointer"
                onClick={() => {
                  if (rule.affected_atoms && rule.affected_atoms.length > 0) {
                    const atomArray = Array.from(molecule.atoms.values())
                    const firstAtomId = atomArray[rule.affected_atoms[0]]?.id
                    if (firstAtomId) {
                      selectAtom(firstAtomId)
                    }
                    if (onHighlightAtoms) {
                      onHighlightAtoms(rule.affected_atoms)
                    }
                  }
                }}
              >
                <div className="flex items-start justify-between mb-2">
                  <div className="flex items-center gap-2">
                    <span className={`px-2 py-1 text-xs font-medium rounded ${getCategoryColor(rule.category)}`}>
                      {rule.category}
                    </span>
                    <span className={`text-xs font-semibold ${getConfidenceColor(rule.confidence)}`}>
                      {(rule.confidence * 100).toFixed(0)}% confidence
                    </span>
                  </div>
                  <button
                    onClick={(e) => {
                      e.stopPropagation()
                      setSelectedRule(rule.rule)
                      setShowExplanation(true)
                    }}
                    className="px-2 py-1 text-xs bg-blue-100 hover:bg-blue-200 text-blue-700 rounded transition-colors"
                  >
                    Explain
                  </button>
                </div>
                <div className="text-sm text-gray-700 mb-1">{rule.description}</div>
                {rule.affected_atoms && rule.affected_atoms.length > 0 && (
                  <div className="text-xs text-gray-500">Affects {rule.affected_atoms.length} atom(s)</div>
                )}
              </div>
            ))}
          </div>
        </div>
      )}

      {/* Export Button */}
      <button
        onClick={() => {
          const dataStr = JSON.stringify(data, null, 2)
          const blob = new Blob([dataStr], { type: 'application/json' })
          const url = URL.createObjectURL(blob)
          const a = document.createElement('a')
          a.href = url
          a.download = 'knowledge_rules.json'
          a.click()
          URL.revokeObjectURL(url)
        }}
        className="w-full px-4 py-2 bg-gray-100 hover:bg-gray-200 rounded transition-colors text-sm"
      >
        Export Knowledge (JSON)
      </button>

      {showExplanation && selectedRule && (
        <PredictionExplanation
          molecule={molecule}
          type="knowledge"
          ruleId={selectedRule}
          onClose={() => {
            setShowExplanation(false)
            setSelectedRule(null)
          }}
        />
      )}
    </div>
  )
}

