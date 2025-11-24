import React from 'react'

interface SharedStatusProps {
  isShared: boolean
  sharedWith?: string[]
  owner?: string
}

export default function SharedStatus({ isShared, sharedWith = [], owner }: SharedStatusProps) {
  if (!isShared) {
    return (
      <div className="p-2 bg-gray-50 rounded text-xs text-gray-600">
        This molecule is private
      </div>
    )
  }

  return (
    <div className="p-2 bg-blue-50 rounded">
      <div className="text-xs font-semibold text-blue-800 mb-1">Shared Molecule</div>
      {owner && (
        <div className="text-xs text-blue-600">Owner: {owner}</div>
      )}
      {sharedWith.length > 0 && (
        <div className="text-xs text-blue-600 mt-1">
          Shared with {sharedWith.length} user(s)
        </div>
      )}
    </div>
  )
}

