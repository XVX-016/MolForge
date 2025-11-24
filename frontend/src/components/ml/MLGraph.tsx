import React, { useState, useEffect } from 'react'
import { MoleculeGraph } from '@biosynth/engine'

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
}

export default function MLGraph({ molecule }: MLGraphProps) {
  const [data, setData] = useState<MLData | null>(null)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)

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
        {data.properties.map((prop, idx) => (
          <div key={idx}>
            <div className="flex items-center justify-between text-xs mb-1">
              <span className="capitalize">{prop.property.replace('_', ' ')}</span>
              <span>{formatValue(prop.value, prop.unit)}</span>
            </div>
            <div className="w-full bg-gray-200 rounded-full h-2">
              <div
                className={`h-2 rounded-full ${
                  prop.property === 'toxicity' ? 'bg-red-500' : 'bg-blue-500'
                }`}
                style={{ width: `${Math.min(100, prop.value * 100)}%` }}
              />
            </div>
          </div>
        ))}
      </div>
    </div>
  )
}

