import React from 'react'

interface AtomContribution {
  atom: number
  contribution: number
  reason: string
  element: string
}

interface ExplanationData {
  property: string
  value: number
  unit: string
  explanation: string
  atom_contributions: AtomContribution[]
}

interface PredictionExplanationProps {
  explanation: ExplanationData
  onClose: () => void
}

export default function PredictionExplanation({
  explanation,
  onClose
}: PredictionExplanationProps) {
  const { property, value, unit, explanation: text, atom_contributions } = explanation

  const formatValue = (val: number, unit: string) => {
    if (unit === 'probability' || unit === 'fraction' || unit === 'score') {
      return `${(val * 100).toFixed(1)}%`
    }
    return `${val.toFixed(2)} ${unit}`
  }

  return (
    <div className="fixed inset-0 bg-black bg-opacity-50 flex items-center justify-center z-50" onClick={onClose}>
      <div
        className="bg-white rounded-lg shadow-xl max-w-2xl w-full mx-4 max-h-[80vh] overflow-y-auto"
        onClick={(e) => e.stopPropagation()}
      >
        <div className="p-6">
          <div className="flex items-center justify-between mb-4">
            <h3 className="text-lg font-semibold text-gray-900">
              Explanation: {property.replace('_', ' ').toUpperCase()}
            </h3>
            <button
              onClick={onClose}
              className="text-gray-400 hover:text-gray-600 transition-colors"
            >
              <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
              </svg>
            </button>
          </div>

          {/* Property Value */}
          <div className="mb-4 p-3 bg-blue-50 rounded-lg">
            <div className="text-sm text-gray-600 mb-1">Predicted Value</div>
            <div className="text-2xl font-bold text-blue-600">{formatValue(value, unit)}</div>
          </div>

          {/* Textual Explanation */}
          <div className="mb-4">
            <h4 className="text-sm font-semibold text-gray-900 mb-2">Explanation</h4>
            <p className="text-gray-700 whitespace-pre-wrap">{text}</p>
          </div>

          {/* Atom Contributions */}
          {atom_contributions && atom_contributions.length > 0 && (
            <div>
              <h4 className="text-sm font-semibold text-gray-900 mb-2">
                Atom Contributions ({atom_contributions.length})
              </h4>
              <div className="space-y-2 max-h-64 overflow-y-auto">
                {atom_contributions.map((contrib, idx) => (
                  <div
                    key={idx}
                    className={`p-3 border rounded-lg ${
                      contrib.contribution > 0
                        ? 'bg-green-50 border-green-200'
                        : 'bg-red-50 border-red-200'
                    }`}
                  >
                    <div className="flex items-center justify-between mb-1">
                      <div className="flex items-center gap-2">
                        <span className="font-semibold text-gray-900">Atom {contrib.atom}</span>
                        <span className="text-xs text-gray-600">({contrib.element})</span>
                      </div>
                      <span
                        className={`font-mono font-semibold ${
                          contrib.contribution > 0 ? 'text-green-700' : 'text-red-700'
                        }`}
                      >
                        {contrib.contribution > 0 ? '+' : ''}{contrib.contribution.toFixed(3)}
                      </span>
                    </div>
                    <div className="text-xs text-gray-600">{contrib.reason}</div>
                  </div>
                ))}
              </div>
            </div>
          )}
        </div>
      </div>
    </div>
  )
}

