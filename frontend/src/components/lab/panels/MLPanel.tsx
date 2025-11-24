import React, { useState } from 'react'
import { useLabStore } from '../../../store/labStore'
import MLGraph from '../../ml/MLGraph'
import MLExport from '../../ml/MLExport'

export default function MLPanel() {
  const molecule = useLabStore(s => s.molecule)
  const [highlightedAtoms, setHighlightedAtoms] = useState<number[]>([])
  const [mlData, setMlData] = useState<any>(null)

  const handleHighlightAtoms = (atomIndices: number[]) => {
    setHighlightedAtoms(atomIndices)
    // TODO: Integrate with MoleculeViewer to highlight atoms
    // This would require passing highlightAtoms prop to MoleculeViewer
  }

  const handleMLDataUpdate = (data: any) => {
    setMlData(data)
  }

  if (!molecule || molecule.atoms.size === 0) {
    return (
      <div className="bg-white rounded-lg border border-gray-200 p-4 shadow-sm">
        <h4 className="text-sm font-semibold text-gray-700 mb-3">ML Predictions</h4>
        <p className="text-xs text-gray-500">Add atoms to run ML predictions</p>
      </div>
    )
  }

  return (
    <div className="bg-white rounded-lg border border-gray-200 p-4 shadow-sm">
      <h4 className="text-sm font-semibold text-gray-700 mb-3">ML Predictions</h4>
      
      <MLGraph
        molecule={molecule}
        onHighlightAtoms={handleHighlightAtoms}
      />
      
      {mlData && (
        <div className="mt-4">
          <MLExport data={mlData} />
        </div>
      )}
      
      {highlightedAtoms.length > 0 && (
        <div className="mt-3 p-2 bg-blue-50 rounded text-xs text-blue-700">
          Highlighting {highlightedAtoms.length} contributing atom(s)
        </div>
      )}
    </div>
  )
}

