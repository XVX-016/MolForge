import React, { useState } from 'react'
import { MoleculeGraph } from '@biosynth/engine'
import { useMoleculeStore } from '../../store/moleculeStore'

interface PathwayStep {
  molecule: any
  step: number
  reaction?: any
  precursors?: any[]
  is_starting?: boolean
}

interface Pathway {
  steps: PathwayStep[]
  total_steps: number
  score?: number
  target: any
}

interface RetrosynthesisStepperProps {
  pathway: Pathway
  onStepSelect?: (step: number) => void
  onMoleculeLoad?: (molecule: any) => void
}

export default function RetrosynthesisStepper({
  pathway,
  onStepSelect,
  onMoleculeLoad
}: RetrosynthesisStepperProps) {
  const [currentStep, setCurrentStep] = useState(0)
  const selectAtom = useMoleculeStore((s) => s.selectAtom)

  if (!pathway || !pathway.steps || pathway.steps.length === 0) {
    return (
      <div className="p-4 text-sm text-gray-500 text-center">
        No pathway data available.
      </div>
    )
  }

  const handleStepClick = (stepIndex: number) => {
    setCurrentStep(stepIndex)
    if (onStepSelect) {
      onStepSelect(stepIndex)
    }
  }

  const handleLoadMolecule = (molecule: any) => {
    if (onMoleculeLoad) {
      onMoleculeLoad(molecule)
    }
  }

  const currentStepData = pathway.steps[currentStep]
  const isFirstStep = currentStep === 0
  const isLastStep = currentStep === pathway.steps.length - 1

  return (
    <div className="space-y-4">
      {/* Step Navigation */}
      <div className="flex items-center justify-between">
        <button
          onClick={() => handleStepClick(Math.max(0, currentStep - 1))}
          disabled={isFirstStep}
          className="px-3 py-1 text-xs bg-gray-100 hover:bg-gray-200 rounded disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
        >
          ← Previous
        </button>
        <div className="text-sm font-semibold text-gray-900">
          Step {currentStep + 1} of {pathway.steps.length}
        </div>
        <button
          onClick={() => handleStepClick(Math.min(pathway.steps.length - 1, currentStep + 1))}
          disabled={isLastStep}
          className="px-3 py-1 text-xs bg-gray-100 hover:bg-gray-200 rounded disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
        >
          Next →
        </button>
      </div>

      {/* Current Step Info */}
      <div className="p-4 bg-blue-50 border border-blue-200 rounded-lg">
        <div className="text-sm font-semibold text-blue-900 mb-2">
          {currentStepData.is_starting ? 'Starting Material' : `Step ${currentStepData.step}`}
        </div>
        {currentStepData.reaction && (
          <div className="text-xs text-blue-700 mb-2">
            Reaction: {currentStepData.reaction.type || 'Unknown'}
          </div>
        )}
        {currentStepData.precursors && currentStepData.precursors.length > 0 && (
          <div className="text-xs text-blue-700">
            Precursors: {currentStepData.precursors.length}
          </div>
        )}
        <button
          onClick={() => handleLoadMolecule(currentStepData.molecule)}
          className="mt-2 px-3 py-1 text-xs bg-blue-600 text-white rounded hover:bg-blue-700 transition-colors"
        >
          Load into Viewer
        </button>
      </div>

      {/* Step Timeline */}
      <div className="space-y-2">
        <div className="text-xs font-semibold text-gray-700 mb-2">Pathway Timeline</div>
        {pathway.steps.map((step, idx) => (
          <div
            key={idx}
            className={`p-2 border rounded cursor-pointer transition-all ${
              idx === currentStep
                ? 'border-blue-500 bg-blue-50'
                : 'border-gray-200 hover:border-gray-300 hover:bg-gray-50'
            }`}
            onClick={() => handleStepClick(idx)}
          >
            <div className="flex items-center justify-between">
              <div className="flex items-center gap-2">
                <div className={`w-2 h-2 rounded-full ${
                  idx === currentStep ? 'bg-blue-600' : 'bg-gray-400'
                }`} />
                <span className="text-xs text-gray-700">
                  {step.is_starting ? 'Starting Material' : `Step ${step.step}`}
                </span>
              </div>
              {step.reaction && (
                <span className="text-xs text-gray-500">
                  {step.reaction.type || 'Reaction'}
                </span>
              )}
            </div>
          </div>
        ))}
      </div>

      {/* Pathway Score */}
      {pathway.score !== undefined && (
        <div className="p-3 bg-gray-50 rounded-lg">
          <div className="text-xs text-gray-600 mb-1">Pathway Score</div>
          <div className="text-lg font-semibold text-gray-900">
            {(pathway.score * 100).toFixed(1)}%
          </div>
          <div className="text-xs text-gray-500 mt-1">
            Total Steps: {pathway.total_steps}
          </div>
        </div>
      )}
    </div>
  )
}

