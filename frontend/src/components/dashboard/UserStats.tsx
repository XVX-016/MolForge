import React from 'react'

interface UserStatsProps {
  data: {
    active_users: number
    total_molecules: number
  }
}

export default function UserStats({ data }: UserStatsProps) {
  const avgMoleculesPerUser = data.active_users > 0 
    ? (data.total_molecules / data.active_users).toFixed(1)
    : '0'

  return (
    <div className="p-4 bg-gradient-to-r from-blue-50 to-purple-50 rounded-lg border border-blue-200">
      <h3 className="text-lg font-semibold text-gray-900 mb-4">User Statistics</h3>
      <div className="grid grid-cols-2 gap-4">
        <div>
          <div className="text-sm text-gray-600 mb-1">Active Users</div>
          <div className="text-3xl font-bold text-blue-600">{data.active_users}</div>
        </div>
        <div>
          <div className="text-sm text-gray-600 mb-1">Avg Molecules/User</div>
          <div className="text-3xl font-bold text-purple-600">{avgMoleculesPerUser}</div>
        </div>
      </div>
    </div>
  )
}

