import React, { useState, useEffect } from 'react'
import { MoleculeGraph } from '@biosynth/engine'
import { useMoleculeStore } from '../../store/moleculeStore'

interface AnalysisData {
  functional_groups: any[]
  molecular_weight: number
  surface_area: number
  complexity: number
  alerts: any[]
  num_atoms: number
  num_bonds: number
  num_rings: number
}

interface AnalysisGraphProps {
  molecule: MoleculeGraph
  onHighlightAtoms?: (atomIndices: number[]) => void
}

export default function AnalysisGraph({ molecule, onHighlightAtoms }: AnalysisGraphProps) {
  const [analysis, setAnalysis] = useState<AnalysisData | null>(null)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const selectAtom = useMoleculeStore((s) => s.selectAtom)

  useEffect(() => {
    if (!molecule || molecule.atoms.size === 0) {
      setAnalysis(null)
      return
    }

    const fetchAnalysis = async () => {
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

        const response = await fetch('/api/kab/analysis', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ atoms, bonds })
        })

        if (!response.ok) {
          throw new Error('Failed to fetch analysis')
        }

        const data = await response.json()
        setAnalysis(data)
      } catch (err) {
        setError(err instanceof Error ? err.message : 'Unknown error')
        console.error('Analysis fetch error:', err)
      } finally {
        setLoading(false)
      }
    }

    fetchAnalysis()
  }, [molecule])

  if (loading) {
    return (
      <div className="flex items-center justify-center py-8">
        <div className="w-6 h-6 border-2 border-blue-600 border-t-transparent rounded-full animate-spin"></div>
        <span className="ml-2 text-sm text-gray-600">Analyzing molecule...</span>
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

  if (!analysis) {
    return (
      <div className="p-4 text-sm text-gray-500 text-center">
        No analysis data available.
      </div>
    )
  }

  const getSeverityColor = (severity: string) => {
    const colors: Record<string, string> = {
      high: 'bg-red-100 text-red-800 border-red-300',
      medium: 'bg-yellow-100 text-yellow-800 border-yellow-300',
      low: 'bg-blue-100 text-blue-800 border-blue-300'
    }
    return colors[severity] || 'bg-gray-100 text-gray-800 border-gray-300'
  }

  return (
    <div className="space-y-6">
      {/* Summary Stats */}
      <div className="grid grid-cols-2 gap-4">
        <div className="p-3 bg-gray-50 rounded-lg">
          <div className="text-xs text-gray-500">Molecular Weight</div>
          <div className="text-lg font-semibold">{analysis.molecular_weight} Da</div>
        </div>
        <div className="p-3 bg-gray-50 rounded-lg">
          <div className="text-xs text-gray-500">Surface Area</div>
          <div className="text-lg font-semibold">{analysis.surface_area} Å²</div>
        </div>
        <div className="p-3 bg-gray-50 rounded-lg">
          <div className="text-xs text-gray-500">Complexity</div>
          <div className="text-lg font-semibold">{analysis.complexity.toFixed(2)}</div>
        </div>
        <div className="p-3 bg-gray-50 rounded-lg">
          <div className="text-xs text-gray-500">Rings</div>
          <div className="text-lg font-semibold">{analysis.num_rings}</div>
        </div>
      </div>

      {/* Functional Groups */}
      {analysis.functional_groups && analysis.functional_groups.length > 0 && (
        <div>
          <h4 className="text-sm font-semibold text-gray-900 mb-2">Functional Groups</h4>
          <div className="space-y-1">
            {analysis.functional_groups.map((group, idx) => (
              <div
                key={idx}
                className="p-2 bg-blue-50 rounded text-sm cursor-pointer hover:bg-blue-100 transition-colors"
                onClick={() => {
                  if (group.atoms && group.atoms.length > 0) {
                    const atomArray = Array.from(molecule.atoms.values())
                    const firstAtomId = atomArray[group.atoms[0]]?.id
                    if (firstAtomId) {
                      selectAtom(firstAtomId)
                    }
                    if (onHighlightAtoms) {
                      onHighlightAtoms(group.atoms)
                    }
                  }
                }}
              >
                <span className="font-medium">{group.type}</span>
                {group.atoms && <span className="text-gray-600 ml-2">({group.atoms.length} atoms)</span>}
              </div>
            ))}
          </div>
        </div>
      )}

      {/* Structural Alerts */}
      {analysis.alerts && analysis.alerts.length > 0 && (
        <div>
          <h4 className="text-sm font-semibold text-gray-900 mb-2">Structural Alerts ({analysis.alerts.length})</h4>
          <div className="space-y-2">
            {analysis.alerts.map((alert, idx) => (
              <div
                key={idx}
                className={`p-3 border rounded-lg ${getSeverityColor(alert.severity)}`}
                onClick={() => {
                  if (onHighlightAtoms && alert.affected_atoms) {
                    onHighlightAtoms(alert.affected_atoms)
                  }
                }}
              >
                <div className="flex items-start justify-between">
                  <div className="flex-1">
                    <div className="font-medium mb-1">{alert.type}</div>
                    <div className="text-sm">{alert.description}</div>
                  </div>
                  <span className={`px-2 py-1 text-xs font-medium rounded ${getSeverityColor(alert.severity)}`}>
                    {alert.severity}
                  </span>
                </div>
              </div>
            ))}
          </div>
        </div>
      )}

      {/* Export Button */}
      <button
        onClick={() => {
          const dataStr = JSON.stringify(analysis, null, 2)
          const blob = new Blob([dataStr], { type: 'application/json' })
          const url = URL.createObjectURL(blob)
          const a = document.createElement('a')
          a.href = url
          a.download = 'molecular_analysis.json'
          a.click()
          URL.revokeObjectURL(url)
        }}
        className="w-full px-4 py-2 bg-gray-100 hover:bg-gray-200 rounded transition-colors text-sm"
      >
        Export Analysis (JSON)
      </button>
    </div>
  )
}

