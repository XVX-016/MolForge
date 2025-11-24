import React from 'react'

interface LabPanelProps {
  children: React.ReactNode
}

/**
 * Centered soft card container for the Lab workspace
 * Matches the target design: beige background, rounded corners, soft shadows
 */
export default function LabPanel({ children }: LabPanelProps) {
  return (
    <div className="w-full flex justify-center mt-8 mb-8">
      <div
        className="rounded-3xl p-6 shadow-[0_4px_20px_rgba(0,0,0,0.08)]"
        style={{
          width: '900px',
          maxWidth: '95vw',
          background: '#f5f1e8', // soft beige
          border: '1px solid #e3dfd6',
        }}
      >
        {children}
      </div>
    </div>
  )
}

