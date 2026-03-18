import React from 'react'

interface PathwayStepReaction {
  name?: string
  type?: string
}

interface PathwayStepMolecule {
  atoms?: unknown[]
}

interface PathwayStepData {
  step?: number
  reaction?: PathwayStepReaction
  molecule?: PathwayStepMolecule
  is_starting?: boolean
}

interface PathwayData {
  score?: number
  total_steps?: number
  steps?: PathwayStepData[]
}

interface PathwayExportProps {
  pathway: PathwayData
  filename?: string
}

export default function PathwayExport({ pathway, filename = 'pathway' }: PathwayExportProps) {
  const handleExport = (format: 'json' | 'csv') => {
    if (format === 'json') {
      const dataStr = JSON.stringify(pathway, null, 2)
      const blob = new Blob([dataStr], { type: 'application/json' })
      const url = URL.createObjectURL(blob)
      const a = document.createElement('a')
      a.href = url
      a.download = `${filename}.json`
      a.click()
      URL.revokeObjectURL(url)
    } else if (format === 'csv') {
      const csv = convertToCSV(pathway)
      const blob = new Blob([csv], { type: 'text/csv' })
      const url = URL.createObjectURL(blob)
      const a = document.createElement('a')
      a.href = url
      a.download = `${filename}.csv`
      a.click()
      URL.revokeObjectURL(url)
    }
  }

  const convertToCSV = (pathway: PathwayData): string => {
    const lines: string[] = []
    
    lines.push('Retrosynthesis Pathway')
    lines.push(`Score,${pathway.score || 'N/A'}`)
    lines.push(`Total Steps,${pathway.total_steps || 0}`)
    lines.push('')
    lines.push('Step,Description,Atoms,Is Starting Material')
    
    if (pathway.steps) {
      pathway.steps.forEach((step, idx) => {
        const desc = step.reaction?.name || step.reaction?.type || `Step ${step.step}`
        const atoms = step.molecule?.atoms?.length || 0
        const isStarting = step.is_starting ? 'Yes' : 'No'
        lines.push(`${step.step || idx},${desc},${atoms},${isStarting}`)
      })
    }
    
    return lines.join('\n')
  }

  return (
    <div className="flex gap-2">
      <button
        onClick={() => handleExport('json')}
        className="px-3 py-1 text-xs bg-blue-100 hover:bg-blue-200 text-blue-700 rounded transition-colors"
      >
        Export JSON
      </button>
      <button
        onClick={() => handleExport('csv')}
        className="px-3 py-1 text-xs bg-green-100 hover:bg-green-200 text-green-700 rounded transition-colors"
      >
        Export CSV
      </button>
    </div>
  )
}
