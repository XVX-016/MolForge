import React from 'react'
import { useLabStore } from '../../../store/labStore'

/**
 * Precision coordinate editing panel
 */
export default function PrecisionPanel() {
  const selectedAtomId = useLabStore(s => s.selectedAtomId)
  const moveAtom = useLabStore(s => s.moveAtom)

  const selectedAtom = selectedAtomId
    ? useLabStore.getState().molecule.atoms.find(a => a.id === selectedAtomId)
    : null

  if (!selectedAtom) {
    return (
      <div className="bg-white rounded-lg border border-gray-200 p-3 shadow-sm">
        <h4 className="text-xs font-semibold text-gray-700 mb-2">Precision Editor</h4>
        <p className="text-xs text-gray-500">Select an atom to edit coordinates</p>
      </div>
    )
  }

  const updatePosition = (axis: 0 | 1 | 2, value: number) => {
    const newPos: [number, number, number] = [...selectedAtom.position]
    newPos[axis] = value
    moveAtom(selectedAtom.id, newPos)
  }

  const nudge = (axis: 0 | 1 | 2, delta: number) => {
    const newPos: [number, number, number] = [...selectedAtom.position]
    newPos[axis] += delta
    moveAtom(selectedAtom.id, newPos)
  }

  return (
    <div className="bg-white rounded-lg border border-gray-200 p-3 shadow-sm">
      <div className="mb-3">
        <h4 className="text-xs font-semibold text-gray-700">Precision Editor</h4>
      </div>

      <div className="space-y-2">
        {(['X', 'Y', 'Z'] as const).map((axis, idx) => (
          <div key={axis} className="flex items-center gap-2">
            <label className="text-xs text-gray-600 w-4">{axis}</label>
            <input
              type="number"
              value={selectedAtom.position[idx].toFixed(3)}
              onChange={(e) => {
                const val = parseFloat(e.target.value) || 0
                updatePosition(idx as 0 | 1 | 2, val)
              }}
              className="flex-1 px-2 py-1 text-xs border border-gray-300 rounded"
              step="0.01"
            />
            <div className="flex gap-1">
              <button
                onClick={() => nudge(idx as 0 | 1 | 2, -0.1)}
                className="px-1.5 py-0.5 text-xs border border-gray-300 rounded hover:bg-gray-50"
                title="-0.1"
              >
                -0.1
              </button>
              <button
                onClick={() => nudge(idx as 0 | 1 | 2, -0.01)}
                className="px-1.5 py-0.5 text-xs border border-gray-300 rounded hover:bg-gray-50"
                title="-0.01"
              >
                -0.01
              </button>
              <button
                onClick={() => nudge(idx as 0 | 1 | 2, 0.01)}
                className="px-1.5 py-0.5 text-xs border border-gray-300 rounded hover:bg-gray-50"
                title="+0.01"
              >
                +0.01
              </button>
              <button
                onClick={() => nudge(idx as 0 | 1 | 2, 0.1)}
                className="px-1.5 py-0.5 text-xs border border-gray-300 rounded hover:bg-gray-50"
                title="+0.1"
              >
                +0.1
              </button>
            </div>
          </div>
        ))}
      </div>
    </div>
  )
}

