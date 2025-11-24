import React, { useState, useEffect } from 'react'
import { MoleculeGraph } from '@biosynth/engine'
import { useMoleculeStore } from '../../store/moleculeStore'
import PredictionExplanation from './PredictionExplanation'

interface MLProperty {
  property: string
  value: number
  unit: string
}

interface MLData {
  properties: MLProperty[]
  model: string
  features?: any
}

interface MLGraphProps {
  molecule: MoleculeGraph
  onHighlightAtoms?: (atomIndices: number[]) => void
}

export default function MLGraph({ molecule, onHighlightAtoms }: MLGraphProps) {
  const [data, setData] = useState<MLData | null>(null)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const [explanations, setExplanations] = useState<Record<string, any>>({})
  const [selectedProperty, setSelectedProperty] = useState<string | null>(null)
  const [showExplanation, setShowExplanation] = useState(false)
  const selectAtom = useMoleculeStore((s) => s.selectAtom)

  useEffect(() => {
    if (!molecule || molecule.atoms.size === 0) {
      setData(null)
      return
    }

    const fetchML = async () => {
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

        const response = await fetch('/api/ml/predict', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ atoms, bonds })
        })

        if (!response.ok) {
          throw new Error('Failed to fetch ML predictions')
        }

        const result = await response.json()
        setData(result)
        
        // Fetch explanations for all properties
        if (result.properties) {
          const explanationPromises = result.properties.map(async (prop: MLProperty) => {
            try {
              const explainRes = await fetch('/api/ml/explain', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                  molecule: { atoms, bonds },
                  property_name: prop.property
                })
              })
              if (explainRes.ok) {
                return await explainRes.json()
              }
            } catch (err) {
              console.error(`Failed to fetch explanation for ${prop.property}:`, err)
            }
            return null
          })
          
          const explanationResults = await Promise.all(explanationPromises)
          const explanationMap: Record<string, any> = {}
          result.properties.forEach((prop: MLProperty, idx: number) => {
            if (explanationResults[idx]) {
              explanationMap[prop.property] = explanationResults[idx]
            }
          })
          setExplanations(explanationMap)
        }
        
        // Notify parent of data update
        if (onHighlightAtoms) {
          // This will be called when user clicks highlight button
        }
      } catch (err) {
        setError(err instanceof Error ? err.message : 'Unknown error')
        console.error('ML fetch error:', err)
      } finally {
        setLoading(false)
      }
    }

    fetchML()
  }, [molecule])

  if (loading) {
    return (
      <div className="flex items-center justify-center py-8">
        <div className="w-6 h-6 border-2 border-blue-600 border-t-transparent rounded-full animate-spin"></div>
        <span className="ml-2 text-sm text-gray-600">Running ML predictions...</span>
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

  if (!data || !data.properties) {
    return (
      <div className="p-4 text-sm text-gray-500 text-center">
        No ML predictions available.
      </div>
    )
  }

  const getPropertyColor = (property: string) => {
    const colors: Record<string, string> = {
      solubility: 'bg-blue-50 text-blue-800',
      toxicity: 'bg-red-50 text-red-800',
      bioavailability: 'bg-green-50 text-green-800',
      drug_likeness: 'bg-purple-50 text-purple-800'
    }
    return colors[property] || 'bg-gray-50 text-gray-800'
  }

  const formatValue = (value: number, unit: string) => {
    if (unit === 'probability' || unit === 'fraction' || unit === 'score') {
      return `${(value * 100).toFixed(1)}%`
    }
    return `${value.toFixed(2)} ${unit}`
  }

  return (
    <div className="space-y-4">
      <div className="flex items-center justify-between mb-4">
        <h4 className="text-sm font-semibold text-gray-900">ML Predictions</h4>
        <span className="text-xs text-gray-500">Model: {data.model}</span>
      </div>

      <div className="grid grid-cols-2 gap-3">
        {data.properties.map((prop, idx) => (
          <div
            key={idx}
            className={`p-3 rounded-lg ${getPropertyColor(prop.property)}`}
          >
            <div className="text-xs font-medium mb-1 capitalize">
              {prop.property.replace('_', ' ')}
            </div>
            <div className="text-lg font-semibold">
              {formatValue(prop.value, prop.unit)}
            </div>
          </div>
        ))}
      </div>

      {/* Progress bars for visual representation */}
      <div className="space-y-2 mt-4">
        {data.properties.map((prop, idx) => {
          const explanation = explanations[prop.property]
          const atomContributions = explanation?.atom_contributions || []
          
          return (
            <div key={idx} className="space-y-2">
              <div className="flex items-center justify-between text-xs mb-1">
                <span className="capitalize">{prop.property.replace('_', ' ')}</span>
                <div className="flex items-center gap-2">
                  <span>{formatValue(prop.value, prop.unit)}</span>
                  {explanation && (
                    <button
                      onClick={() => {
                        setSelectedProperty(prop.property)
                        setShowExplanation(true)
                      }}
                      className="px-2 py-0.5 text-xs bg-blue-100 hover:bg-blue-200 text-blue-700 rounded transition-colors"
                    >
                      Explain
                    </button>
                  )}
                </div>
              </div>
              <div className="w-full bg-gray-200 rounded-full h-2">
                <div
                  className={`h-2 rounded-full ${
                    prop.property === 'toxicity' ? 'bg-red-500' : 'bg-blue-500'
                  }`}
                  style={{ width: `${Math.min(100, prop.value * 100)}%` }}
                />
              </div>
              
              {/* Atom Contributions */}
              {atomContributions.length > 0 && (
                <div className="mt-2">
                  <button
                    onClick={() => {
                      const contributingAtoms = atomContributions
                        .slice(0, 5)
                        .map((c: any) => c.atom)
                      const atomArray = Array.from(molecule.atoms.values())
                      const firstAtomId = atomArray[contributingAtoms[0]]?.id
                      if (firstAtomId) {
                        selectAtom(firstAtomId)
                      }
                      if (onHighlightAtoms) {
                        onHighlightAtoms(contributingAtoms)
                      }
                    }}
                    className="text-xs text-purple-600 hover:text-purple-800"
                  >
                    Highlight contributing atoms ({atomContributions.length})
                  </button>
                  <div className="mt-1 space-y-1 max-h-20 overflow-y-auto">
                    {atomContributions.slice(0, 3).map((contrib: any, cIdx: number) => (
                      <div key={cIdx} className="text-xs text-gray-600">
                        Atom {contrib.atom} ({contrib.element}): {contrib.contribution > 0 ? '+' : ''}{contrib.contribution.toFixed(3)} - {contrib.reason}
                      </div>
                    ))}
                  </div>
                </div>
              )}
            </div>
          )
        })}
      </div>
      
      {showExplanation && selectedProperty && explanations[selectedProperty] && (
        <PredictionExplanation
          explanation={explanations[selectedProperty]}
          onClose={() => {
            setShowExplanation(false)
            setSelectedProperty(null)
          }}
        />
      )}
    </div>
  )
}

