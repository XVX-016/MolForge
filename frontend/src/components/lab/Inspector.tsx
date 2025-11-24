import React from 'react'
import { useLabStore } from '../../store/labStore'
import { VALENCE_LIMITS } from '../../utils/bondRules'

export default function Inspector() {
  const mol = useLabStore(s => s.molecule)
  const selectedAtomId = useLabStore(s => s.selectedAtomId)
  const auto = useLabStore(s => s.autoBond)
  const setAuto = useLabStore(s => s.setAutoBond)
  const undo = useLabStore(s => s.undo)
  const redo = useLabStore(s => s.redo)
  const resetMolecule = useLabStore(s => s.resetMolecule)
  const moveAtom = useLabStore(s => s.moveAtom)
  
  const selectedAtom = selectedAtomId 
    ? mol.atoms.find(a => a.id === selectedAtomId)
    : null

  // Count bonds for selected atom
  const bondCount = selectedAtom
    ? mol.bonds.filter(b => b.atom1 === selectedAtom.id || b.atom2 === selectedAtom.id).length
    : 0

  const maxValence = selectedAtom ? (VALENCE_LIMITS[selectedAtom.element] || 4) : 0
  
  return (
    <div className="bg-white rounded-lg border border-gray-200 p-3 shadow-sm">
      <h4 className="text-xs font-semibold text-gray-700 mb-3 text-center">Inspector</h4>
      <div className="space-y-2 text-xs">
        <div className="flex justify-between">
          <span className="text-gray-600">Atoms:</span>
          <span className="font-semibold text-gray-900">{mol.atoms.length}</span>
        </div>
        <div className="flex justify-between">
          <span className="text-gray-600">Bonds:</span>
          <span className="font-semibold text-gray-900">{mol.bonds.length}</span>
        </div>
      </div>

      {/* Selected Atom Properties */}
      {selectedAtom && (
        <div className="mt-3 pt-3 border-t border-gray-200 space-y-2">
          <div className="text-xs font-semibold text-gray-700 mb-2">Selected Atom</div>
          <div className="space-y-1.5 text-xs">
            <div className="flex justify-between">
              <span className="text-gray-600">Element:</span>
              <span className="font-semibold text-gray-900">{selectedAtom.element}</span>
            </div>
            <div className="flex justify-between">
              <span className="text-gray-600">Bonds:</span>
              <span className="font-semibold text-gray-900">
                {bondCount} / {maxValence}
              </span>
            </div>
            <div className="flex justify-between">
              <span className="text-gray-600">Charge:</span>
              <input
                type="number"
                value={selectedAtom.charge || 0}
                onChange={(e) => {
                  const charge = parseInt(e.target.value) || 0
                  const updatedAtom = { ...selectedAtom, charge }
                  const updatedAtoms = mol.atoms.map(a => 
                    a.id === selectedAtom.id ? updatedAtom : a
                  )
                  useLabStore.setState({ 
                    molecule: { ...mol, atoms: updatedAtoms } 
                  })
                }}
                className="w-12 px-1 py-0.5 text-xs border border-gray-300 rounded text-right"
                min="-3"
                max="3"
              />
            </div>
            <div className="text-[10px] text-gray-500 mt-2">
              Position: ({selectedAtom.position[0].toFixed(2)}, {selectedAtom.position[1].toFixed(2)}, {selectedAtom.position[2].toFixed(2)})
            </div>
          </div>
        </div>
      )}
      <div className="mt-3 pt-3 border-t border-gray-200">
        <label className="flex items-center gap-2 text-xs cursor-pointer">
          <input
            type="checkbox"
            checked={auto}
            onChange={(e) => setAuto(e.target.checked)}
            className="rounded"
          />
          <span className="text-gray-700">Auto-bond</span>
        </label>
      </div>
      <div className="mt-3 pt-3 border-t border-gray-200 grid grid-cols-3 gap-1">
        <button
          onClick={undo}
          className="px-2 py-1 text-[10px] rounded border border-gray-300 bg-white hover:bg-gray-50 transition-colors"
          title="Undo"
        >
          ↶
        </button>
        <button
          onClick={redo}
          className="px-2 py-1 text-[10px] rounded border border-gray-300 bg-white hover:bg-gray-50 transition-colors"
          title="Redo"
        >
          ↷
        </button>
        <button
          onClick={resetMolecule}
          className="px-2 py-1 text-[10px] rounded border border-red-300 bg-white hover:bg-red-50 text-red-600 transition-colors"
          title="Reset"
        >
          ×
        </button>
      </div>
    </div>
  )
}

