import React from 'react'

interface DashboardChartsProps {
  data: {
    total_molecules: number
    molecules_with_predictions: number
    molecules_with_3d: number
    spectra_computed: number
    energy_calculations: number
    kab_analyses: number
    active_users: number
  }
}

export default function DashboardCharts({ data }: DashboardChartsProps) {
  return (
    <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
      {/* Total Molecules */}
      <div className="p-4 bg-white border border-gray-200 rounded-lg shadow-sm">
        <div className="text-sm text-gray-600 mb-1">Total Molecules</div>
        <div className="text-2xl font-bold text-gray-900">{data.total_molecules}</div>
      </div>

      {/* Molecules with Predictions */}
      <div className="p-4 bg-white border border-gray-200 rounded-lg shadow-sm">
        <div className="text-sm text-gray-600 mb-1">With Predictions</div>
        <div className="text-2xl font-bold text-blue-600">{data.molecules_with_predictions}</div>
      </div>

      {/* Molecules with 3D */}
      <div className="p-4 bg-white border border-gray-200 rounded-lg shadow-sm">
        <div className="text-sm text-gray-600 mb-1">3D Structures</div>
        <div className="text-2xl font-bold text-green-600">{data.molecules_with_3d}</div>
      </div>

      {/* Spectra Computed */}
      <div className="p-4 bg-white border border-gray-200 rounded-lg shadow-sm">
        <div className="text-sm text-gray-600 mb-1">Spectra Computed</div>
        <div className="text-2xl font-bold text-purple-600">{data.spectra_computed}</div>
      </div>

      {/* Energy Calculations */}
      <div className="p-4 bg-white border border-gray-200 rounded-lg shadow-sm">
        <div className="text-sm text-gray-600 mb-1">Energy Calculations</div>
        <div className="text-2xl font-bold text-orange-600">{data.energy_calculations}</div>
      </div>

      {/* KAB Analyses */}
      <div className="p-4 bg-white border border-gray-200 rounded-lg shadow-sm">
        <div className="text-sm text-gray-600 mb-1">KAB Analyses</div>
        <div className="text-2xl font-bold text-red-600">{data.kab_analyses}</div>
      </div>
    </div>
  )
}

