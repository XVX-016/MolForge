import React from 'react'

interface CollabExportProps {
  moleculeId: string
  onExport?: () => void
}

export default function CollabExport({ moleculeId, onExport }: CollabExportProps) {
  const handleExport = async (format: 'json' | 'mol') => {
    try {
      const response = await fetch(`/api/collaboration/load/${moleculeId}`)
      if (!response.ok) {
        throw new Error('Failed to load molecule')
      }

      const data = await response.json()
      
      if (format === 'json') {
        const dataStr = JSON.stringify(data, null, 2)
        const blob = new Blob([dataStr], { type: 'application/json' })
        const url = URL.createObjectURL(blob)
        const a = document.createElement('a')
        a.href = url
        a.download = `molecule_${moleculeId}.json`
        a.click()
        URL.revokeObjectURL(url)
      }

      if (onExport) {
        onExport()
      }
    } catch (err) {
      alert(err instanceof Error ? err.message : 'Export failed')
    }
  }

  return (
    <div className="flex gap-2">
      <button
        onClick={() => handleExport('json')}
        className="px-3 py-1 text-xs bg-blue-100 hover:bg-blue-200 text-blue-700 rounded transition-colors"
      >
        Export JSON
      </button>
    </div>
  )
}

