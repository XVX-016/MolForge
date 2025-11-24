import React, { useState, useEffect } from 'react'
import DashboardCharts from './DashboardCharts'
import DashboardFilters from './DashboardFilters'
import UserStats from './UserStats'

interface DashboardData {
  total_molecules: number
  molecules_with_predictions: number
  molecules_with_3d: number
  spectra_computed: number
  energy_calculations: number
  kab_analyses: number
  active_users: number
}

export default function Dashboard() {
  const [data, setData] = useState<DashboardData | null>(null)
  const [loading, setLoading] = useState(true)
  const [filters, setFilters] = useState({})

  useEffect(() => {
    const fetchDashboard = async () => {
      try {
        const response = await fetch('/api/dashboard/summary')
        if (!response.ok) {
          throw new Error('Failed to fetch dashboard data')
        }
        const result = await response.json()
        setData(result)
      } catch (err) {
        console.error('Dashboard fetch error:', err)
      } finally {
        setLoading(false)
      }
    }

    fetchDashboard()
  }, [])

  if (loading) {
    return (
      <div className="flex items-center justify-center py-12">
        <div className="w-8 h-8 border-2 border-blue-600 border-t-transparent rounded-full animate-spin"></div>
      </div>
    )
  }

  if (!data) {
    return (
      <div className="p-4 text-sm text-gray-500 text-center">
        Failed to load dashboard data
      </div>
    )
  }

  return (
    <div className="space-y-6">
      <div className="flex items-center justify-between">
        <h2 className="text-2xl font-bold text-gray-900">Dashboard</h2>
        <DashboardFilters filters={filters} onFiltersChange={setFilters} />
      </div>

      <UserStats data={data} />

      <DashboardCharts data={data} />
    </div>
  )
}

