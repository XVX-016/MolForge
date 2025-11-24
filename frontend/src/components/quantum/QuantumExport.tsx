import React from 'react'

interface QuantumExportProps {
  data: any
  filename?: string
}

export default function QuantumExport({ data, filename = 'quantum_data' }: QuantumExportProps) {
  const handleExport = (format: 'json' | 'csv') => {
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
      // Convert to CSV format
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
    // Simple CSV conversion
    const lines: string[] = []
    if (data.homo_lumo) {
      lines.push('Property,Value,Unit')
      lines.push(`HOMO,${data.homo_lumo.HOMO},eV`)
      lines.push(`LUMO,${data.homo_lumo.LUMO},eV`)
      lines.push(`Gap,${data.homo_lumo.gap},eV`)
    }
    if (data.esp && data.esp.esp_values) {
      lines.push('')
      lines.push('Atom,ESP')
      data.esp.esp_values.forEach((item: any) => {
        lines.push(`${item.atom},${item.esp}`)
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

