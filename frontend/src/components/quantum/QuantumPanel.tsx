import React, { useState, useEffect } from 'react'
import { MoleculeGraph } from '@biosynth/engine'
import { useMoleculeStore } from '../../store/moleculeStore'

interface QuantumData {
  homo_lumo: {
    HOMO: number
    LUMO: number
    gap: number
    contributing_atoms: number[]
  }
  esp: {
    esp_values: Array<{ atom: number; esp: number }>
  }
}

interface QuantumPanelProps {
  molecule: MoleculeGraph
  onHighlightAtoms?: (atomIndices: number[]) => void
}

export default function QuantumPanel({ molecule, onHighlightAtoms }: QuantumPanelProps) {
  const [data, setData] = useState<QuantumData | null>(null)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const selectAtom = useMoleculeStore((s) => s.selectAtom)

  useEffect(() => {
    if (!molecule || molecule.atoms.size === 0) {
      setData(null)
      return
    }

    const fetchQuantum = async () => {
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

        const response = await fetch('/api/quantum/calculate', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ atoms, bonds })
        })

        if (!response.ok) {
          throw new Error('Failed to fetch quantum data')
        }

        const result = await response.json()
        setData(result)
      } catch (err) {
        setError(err instanceof Error ? err.message : 'Unknown error')
        console.error('Quantum fetch error:', err)
      } finally {
        setLoading(false)
      }
    }

    fetchQuantum()
  }, [molecule])

  if (loading) {
    return (
      <div className="flex items-center justify-center py-8">
        <div className="w-6 h-6 border-2 border-blue-600 border-t-transparent rounded-full animate-spin"></div>
        <span className="ml-2 text-sm text-gray-600">Calculating quantum properties...</span>
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
        No quantum data available.
      </div>
    )
  }

  const { homo_lumo, esp } = data

  return (
    <div className="space-y-6">
      {/* HOMO/LUMO */}
      <div>
        <h4 className="text-sm font-semibold text-gray-900 mb-3">Orbital Energies</h4>
        <div className="grid grid-cols-3 gap-3 mb-3">
          <div className="p-3 bg-blue-50 rounded-lg">
            <div className="text-xs text-gray-600 mb-1">HOMO</div>
            <div className="text-lg font-semibold">{homo_lumo.HOMO} eV</div>
            <div className="text-xs text-gray-500 mt-1">Highest Occupied</div>
          </div>
          <div className="p-3 bg-green-50 rounded-lg">
            <div className="text-xs text-gray-600 mb-1">LUMO</div>
            <div className="text-lg font-semibold">{homo_lumo.LUMO} eV</div>
            <div className="text-xs text-gray-500 mt-1">Lowest Unoccupied</div>
          </div>
          <div className="p-3 bg-purple-50 rounded-lg">
            <div className="text-xs text-gray-600 mb-1">Gap</div>
            <div className="text-lg font-semibold">{homo_lumo.gap} eV</div>
            <div className="text-xs text-gray-500 mt-1">Energy Gap</div>
          </div>
        </div>
        
        {/* Orbital Visualization */}
        <div className="p-3 bg-gray-50 rounded-lg mb-3">
          <div className="text-xs text-gray-600 mb-2">Orbital Energy Diagram</div>
          <div className="relative h-24 border-l-2 border-gray-400">
            {/* HOMO line */}
            <div 
              className="absolute left-0 right-0 border-t-2 border-blue-500"
              style={{ top: '30%' }}
            >
              <div className="absolute -left-8 text-xs text-blue-600 font-semibold">HOMO</div>
              <div className="absolute -right-12 text-xs text-blue-600">{homo_lumo.HOMO} eV</div>
            </div>
            {/* LUMO line */}
            <div 
              className="absolute left-0 right-0 border-t-2 border-green-500"
              style={{ top: '70%' }}
            >
              <div className="absolute -left-8 text-xs text-green-600 font-semibold">LUMO</div>
              <div className="absolute -right-12 text-xs text-green-600">{homo_lumo.LUMO} eV</div>
            </div>
            {/* Gap arrow */}
            <div 
              className="absolute left-1/2 transform -translate-x-1/2"
              style={{ top: '30%', height: '40%' }}
            >
              <div className="h-full border-l-2 border-dashed border-purple-500"></div>
              <div className="absolute -bottom-4 left-1/2 transform -translate-x-1/2 text-xs text-purple-600 font-semibold">
                {homo_lumo.gap} eV gap
              </div>
            </div>
          </div>
        </div>
        
        {homo_lumo.contributing_atoms && homo_lumo.contributing_atoms.length > 0 && (
          <button
            onClick={() => {
              const atomArray = Array.from(molecule.atoms.values())
              const firstAtomId = atomArray[homo_lumo.contributing_atoms[0]]?.id
              if (firstAtomId) {
                selectAtom(firstAtomId)
              }
              if (onHighlightAtoms) {
                onHighlightAtoms(homo_lumo.contributing_atoms)
              }
            }}
            className="w-full px-3 py-2 text-xs bg-blue-100 hover:bg-blue-200 text-blue-700 rounded transition-colors"
          >
            Highlight π System ({homo_lumo.contributing_atoms.length} atoms)
          </button>
        )}
      </div>

      {/* ESP */}
      {esp.esp_values && esp.esp_values.length > 0 && (
        <div>
          <h4 className="text-sm font-semibold text-gray-900 mb-3">Electrostatic Potential</h4>
          
          {/* ESP Color Scale */}
          <div className="mb-3 p-2 bg-gray-50 rounded-lg">
            <div className="text-xs text-gray-600 mb-2">ESP Scale</div>
            <div className="flex items-center gap-2">
              <div className="flex-1 h-4 bg-gradient-to-r from-red-500 via-yellow-500 to-blue-500 rounded"></div>
              <div className="text-xs text-gray-500">Negative → Positive</div>
            </div>
          </div>
          
          {/* ESP Values List */}
          <div className="space-y-2 max-h-48 overflow-y-auto">
            {esp.esp_values
              .sort((a, b) => b.esp - a.esp) // Sort by ESP value
              .slice(0, 10)
              .map((item, idx) => {
                const espValue = item.esp
                const isNegative = espValue < 0
                const intensity = Math.min(1, Math.abs(espValue) / 0.5) // Normalize to 0-1
                
                return (
                  <div
                    key={idx}
                    className="flex items-center justify-between p-2 bg-gray-50 rounded text-sm cursor-pointer hover:bg-gray-100 transition-colors"
                    onClick={() => {
                      const atomArray = Array.from(molecule.atoms.values())
                      const atomId = atomArray[item.atom]?.id
                      if (atomId) {
                        selectAtom(atomId)
                      }
                      if (onHighlightAtoms) {
                        onHighlightAtoms([item.atom])
                      }
                    }}
                  >
                    <div className="flex items-center gap-2">
                      <div
                        className={`w-3 h-3 rounded-full ${
                          isNegative 
                            ? `bg-red-${Math.min(500, 300 + Math.floor(intensity * 200))}` 
                            : `bg-blue-${Math.min(500, 300 + Math.floor(intensity * 200))}`
                        }`}
                        style={{
                          backgroundColor: isNegative
                            ? `rgba(239, 68, 68, ${0.5 + intensity * 0.5})`
                            : `rgba(59, 130, 246, ${0.5 + intensity * 0.5})`
                        }}
                      />
                      <span>Atom {item.atom}</span>
                    </div>
                    <span className={`font-mono font-semibold ${
                      isNegative ? 'text-red-600' : 'text-blue-600'
                    }`}>
                      {item.esp > 0 ? '+' : ''}{item.esp.toFixed(3)}
                    </span>
                  </div>
                )
              })}
            {esp.esp_values.length > 10 && (
              <div className="text-xs text-gray-500 text-center">
                +{esp.esp_values.length - 10} more atoms
              </div>
            )}
          </div>
          
          <button
            onClick={() => {
              if (onHighlightAtoms) {
                // Highlight atoms with extreme ESP values
                const sorted = [...esp.esp_values].sort((a, b) => Math.abs(b.esp) - Math.abs(a.esp))
                const extremeAtoms = sorted.slice(0, 5).map(item => item.atom)
                onHighlightAtoms(extremeAtoms)
              }
            }}
            className="mt-2 w-full px-3 py-2 text-xs bg-purple-100 hover:bg-purple-200 text-purple-700 rounded transition-colors"
          >
            Highlight Extreme ESP Values
          </button>
        </div>
      )}

      {/* Export */}
      <button
        onClick={() => {
          const dataStr = JSON.stringify(data, null, 2)
          const blob = new Blob([dataStr], { type: 'application/json' })
          const url = URL.createObjectURL(blob)
          const a = document.createElement('a')
          a.href = url
          a.download = 'quantum_properties.json'
          a.click()
          URL.revokeObjectURL(url)
        }}
        className="w-full px-4 py-2 bg-gray-100 hover:bg-gray-200 rounded transition-colors text-sm"
      >
        Export Quantum Data (JSON)
      </button>
    </div>
  )
}

