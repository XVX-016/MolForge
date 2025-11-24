import React from 'react'

interface ReactionExportProps {
  reactionData: any
  mechanismData?: any
  filename?: string
}

export default function ReactionExport({
  reactionData,
  mechanismData,
  filename = 'reaction_data'
}: ReactionExportProps) {
  const handleExport = (format: 'json' | 'csv') => {
    const data = {
      reaction: reactionData,
      mechanism: mechanismData
    }

    if (format === 'json') {
      const dataStr = JSON.stringify(data, null, 2)
      const blob = new Blob([dataStr], { type: 'application/json' })
      const url = URL.createObjectURL(blob)
      const a = document.createElement('a')
      a.href = url
      a.download = `${filename}.json`
      a.click()
      URL.revokeObjectURL(url)
    } else if (format === 'csv') {
      const csv = convertToCSV(data)
      const blob = new Blob([csv], { type: 'text/csv' })
      const url = URL.createObjectURL(blob)
      const a = document.createElement('a')
      a.href = url
      a.download = `${filename}.csv`
      a.click()
      URL.revokeObjectURL(url)
    }
  }

  const convertToCSV = (data: any): string => {
    const lines: string[] = []
    
    if (data.reaction) {
      lines.push('Reaction Data')
      lines.push('Property,Value')
      if (data.reaction.reaction_type) {
        lines.push(`Reaction Type,${data.reaction.reaction_type}`)
      }
      if (data.reaction.products) {
        lines.push(`Number of Products,${data.reaction.products.length}`)
      }
      lines.push('')
    }

    if (data.mechanism && data.mechanism.steps) {
      lines.push('Mechanism Steps')
      lines.push('Step,Description,Energy (kcal/mol)')
      data.mechanism.steps.forEach((step: any) => {
        lines.push(`${step.step},${step.description},${step.energy}`)
      })
      lines.push('')
      if (data.mechanism.activation_energy !== undefined) {
        lines.push(`Activation Energy,${data.mechanism.activation_energy}`)
      }
      if (data.mechanism.reaction_energy !== undefined) {
        lines.push(`Reaction Energy,${data.mechanism.reaction_energy}`)
      }
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

