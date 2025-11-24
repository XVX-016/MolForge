import React, { useState } from 'react'
import { MoleculeGraph } from '@biosynth/engine'
import { moleculeToJSON } from '../../lib/engineAdapter'

interface CollabPanelProps {
  molecule: MoleculeGraph
  userId?: string
}

export default function CollabPanel({ molecule, userId }: CollabPanelProps) {
  const [saving, setSaving] = useState(false)
  const [moleculeName, setMoleculeName] = useState('')
  const [saved, setSaved] = useState(false)

  const handleSave = async () => {
    if (!molecule || !userId) {
      alert('Please sign in to save molecules')
      return
    }

    if (!moleculeName.trim()) {
      alert('Please enter a molecule name')
      return
    }

    setSaving(true)
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

      const response = await fetch('/api/collaboration/save', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          molecule: { atoms, bonds },
          user_id: userId,
          molecule_name: moleculeName
        })
      })

      if (!response.ok) {
        throw new Error('Failed to save molecule')
      }

      const result = await response.json()
      setSaved(true)
      setTimeout(() => setSaved(false), 3000)
    } catch (err) {
      alert(err instanceof Error ? err.message : 'Failed to save molecule')
    } finally {
      setSaving(false)
    }
  }

  if (!userId) {
    return (
      <div className="p-4 text-sm text-gray-500 text-center">
        Sign in to use collaboration features
      </div>
    )
  }

  return (
    <div className="space-y-4">
      <h4 className="text-sm font-semibold text-gray-900">Cloud Storage</h4>
      
      <div>
        <label className="block text-xs text-gray-600 mb-1">Molecule Name</label>
        <input
          type="text"
          value={moleculeName}
          onChange={(e) => setMoleculeName(e.target.value)}
          placeholder="Enter molecule name"
          className="w-full px-3 py-2 text-sm border border-gray-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500"
        />
      </div>

      <button
        onClick={handleSave}
        disabled={saving || !moleculeName.trim()}
        className="w-full px-4 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700 disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
      >
        {saving ? 'Saving...' : saved ? 'Saved!' : 'Save to Cloud'}
      </button>

      <div className="pt-4 border-t border-gray-200">
        <h5 className="text-xs font-semibold text-gray-700 mb-2">Features</h5>
        <ul className="text-xs text-gray-600 space-y-1">
          <li>• Save molecules to cloud</li>
          <li>• Share with team members</li>
          <li>• Version history tracking</li>
          <li>• Fork and collaborate</li>
        </ul>
      </div>
    </div>
  )
}

