import React, { useState, useEffect } from 'react'
import { MoleculeGraph } from '@biosynth/engine'
import { moleculeToJSON } from '../../lib/engineAdapter'
import { useMoleculeStore } from '../../store/moleculeStore'
import PredictionExplanation from './PredictionExplanation'

interface BindingSite {
  atom_indices: number[]
  score: number
  type: string
  center: [number, number, number]
  group_type?: string
}

interface BindingGraphProps {
  molecule: MoleculeGraph
  onHighlightAtoms?: (atomIndices: number[]) => void
}

export default function BindingGraph({ molecule, onHighlightAtoms }: BindingGraphProps) {
  const [sites, setSites] = useState<BindingSite[]>([])
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const [selectedSite, setSelectedSite] = useState<number | null>(null)
  const [showExplanation, setShowExplanation] = useState(false)
  const selectAtom = useMoleculeStore((s) => s.selectAtom)

  useEffect(() => {
    if (!molecule || molecule.atoms.size === 0) {
      setSites([])
      return
    }

    const fetchBindingSites = async () => {
      setLoading(true)
      setError(null)
      try {
        // Convert molecule to API format
        const jsonGraph = moleculeToJSON(molecule)
        const molData = JSON.parse(jsonGraph)
        
        // Convert to API format (atoms and bonds arrays)
        const atoms = Array.from(molecule.atoms.values()).map((atom, idx) => ({
          element: atom.element,
          position: atom.position,
          id: idx
        }))
        
        const bonds: any[] = []
        molecule.bonds.forEach((bond, idx) => {
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

        const response = await fetch('/api/kab/binding', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ atoms, bonds })
        })

        if (!response.ok) {
          throw new Error('Failed to fetch binding sites')
        }

        const data = await response.json()
        setSites(data.sites || [])
      } catch (err) {
        setError(err instanceof Error ? err.message : 'Unknown error')
        console.error('Binding sites fetch error:', err)
      } finally {
        setLoading(false)
      }
    }

    fetchBindingSites()
  }, [molecule])

  const handleSiteClick = (siteIndex: number) => {
    setSelectedSite(siteIndex)
    const site = sites[siteIndex]
    if (site && site.atom_indices.length > 0) {
      // Map atom indices to atom IDs
      const atomArray = Array.from(molecule.atoms.values())
      const firstAtomId = atomArray[site.atom_indices[0]]?.id
      if (firstAtomId) {
        selectAtom(firstAtomId)
      }
      if (onHighlightAtoms) {
        onHighlightAtoms(site.atom_indices)
      }
    }
  }

  const handleExplain = (siteIndex: number) => {
    setSelectedSite(siteIndex)
    setShowExplanation(true)
  }

  const getScoreColor = (score: number) => {
    if (score >= 0.7) return 'text-green-600 bg-green-50'
    if (score >= 0.4) return 'text-yellow-600 bg-yellow-50'
    return 'text-red-600 bg-red-50'
  }

  const getTypeColor = (type: string) => {
    const colors: Record<string, string> = {
      polar: 'bg-blue-100 text-blue-800',
      hydrophobic: 'bg-gray-100 text-gray-800',
      aromatic: 'bg-purple-100 text-purple-800'
    }
    return colors[type] || 'bg-gray-100 text-gray-800'
  }

  if (loading) {
    return (
      <div className="flex items-center justify-center py-8">
        <div className="w-6 h-6 border-2 border-blue-600 border-t-transparent rounded-full animate-spin"></div>
        <span className="ml-2 text-sm text-gray-600">Analyzing binding sites...</span>
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

  if (sites.length === 0) {
    return (
      <div className="p-4 text-sm text-gray-500 text-center">
        No binding sites detected. Add functional groups to see predictions.
      </div>
    )
  }

  return (
    <div className="space-y-4">
      <div className="flex items-center justify-between mb-4">
        <h3 className="text-lg font-semibold text-gray-900">
          Binding Sites ({sites.length})
        </h3>
        <button
          onClick={() => {
            const dataStr = JSON.stringify(sites, null, 2)
            const blob = new Blob([dataStr], { type: 'application/json' })
            const url = URL.createObjectURL(blob)
            const a = document.createElement('a')
            a.href = url
            a.download = 'binding_sites.json'
            a.click()
            URL.revokeObjectURL(url)
          }}
          className="px-3 py-1 text-xs bg-gray-100 hover:bg-gray-200 rounded transition-colors"
        >
          Export JSON
        </button>
      </div>

      <div className="space-y-2">
        {sites.map((site, index) => (
          <div
            key={index}
            className={`p-3 border rounded-lg cursor-pointer transition-all ${
              selectedSite === index
                ? 'border-blue-500 bg-blue-50'
                : 'border-gray-200 hover:border-gray-300 hover:bg-gray-50'
            }`}
            onClick={() => handleSiteClick(index)}
          >
            <div className="flex items-start justify-between">
              <div className="flex-1">
                <div className="flex items-center gap-2 mb-2">
                  <span className={`px-2 py-1 text-xs font-medium rounded ${getTypeColor(site.type)}`}>
                    {site.type}
                  </span>
                  <span className={`px-2 py-1 text-xs font-semibold rounded ${getScoreColor(site.score)}`}>
                    Score: {(site.score * 100).toFixed(1)}%
                  </span>
                </div>
                <div className="text-sm text-gray-600">
                  <p>Atoms: {site.atom_indices.length}</p>
                  {site.group_type && (
                    <p className="text-xs text-gray-500 mt-1">Group: {site.group_type}</p>
                  )}
                </div>
              </div>
              <button
                onClick={(e) => {
                  e.stopPropagation()
                  handleExplain(index)
                }}
                className="ml-2 px-2 py-1 text-xs bg-blue-100 hover:bg-blue-200 text-blue-700 rounded transition-colors"
              >
                Explain
              </button>
            </div>
          </div>
        ))}
      </div>

      {showExplanation && selectedSite !== null && (
        <PredictionExplanation
          molecule={molecule}
          type="binding"
          siteId={selectedSite}
          onClose={() => {
            setShowExplanation(false)
            setSelectedSite(null)
          }}
        />
      )}
    </div>
  )
}

