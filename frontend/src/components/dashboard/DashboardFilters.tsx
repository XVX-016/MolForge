import React from 'react'

interface DashboardFiltersProps {
  filters: Record<string, any>
  onFiltersChange: (filters: Record<string, any>) => void
}

export default function DashboardFilters({ filters, onFiltersChange }: DashboardFiltersProps) {
  return (
    <div className="flex gap-2">
      <select
        className="px-3 py-1 text-sm border border-gray-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500"
        onChange={(e) => onFiltersChange({ ...filters, timeRange: e.target.value })}
      >
        <option value="all">All Time</option>
        <option value="week">Last Week</option>
        <option value="month">Last Month</option>
        <option value="year">Last Year</option>
      </select>
    </div>
  )
}

