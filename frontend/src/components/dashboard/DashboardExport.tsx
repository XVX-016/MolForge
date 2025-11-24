import React from 'react'

interface DashboardExportProps {
  data: any
}

export default function DashboardExport({ data }: DashboardExportProps) {
  const handleExport = (format: 'json' | 'csv') => {
    if (format === 'json') {
      const dataStr = JSON.stringify(data, null, 2)
      const blob = new Blob([dataStr], { type: 'application/json' })
      const url = URL.createObjectURL(blob)
      const a = document.createElement('a')
      a.href = url
      a.download = 'dashboard_data.json'
      a.click()
      URL.revokeObjectURL(url)
    } else if (format === 'csv') {
      const csv = convertToCSV(data)
      const blob = new Blob([csv], { type: 'text/csv' })
      const url = URL.createObjectURL(blob)
      const a = document.createElement('a')
      a.href = url
      a.download = 'dashboard_data.csv'
      a.click()
      URL.revokeObjectURL(url)
    }
  }

  const convertToCSV = (data: any): string => {
    const lines: string[] = ['Metric,Value']
    Object.entries(data).forEach(([key, value]) => {
      lines.push(`${key},${value}`)
    })
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

