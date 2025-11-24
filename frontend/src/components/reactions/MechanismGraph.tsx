import React, { useState } from 'react'
import { MoleculeGraph } from '@biosynth/engine'
import { useMoleculeStore } from '../../store/moleculeStore'
import { useLabStore } from '../../store/labStore'

interface MechanismStep {
  step: number
  intermediate: any
  energy: number
  description: string
}

interface MechanismData {
  steps: MechanismStep[]
  activation_energy: number
  reaction_energy: number
}

interface MechanismGraphProps {
  molecule?: MoleculeGraph | any
  mechanismData: MechanismData | null
  onStepSelect?: (step: number) => void
  onHighlightAtoms?: (atomIndices: number[]) => void
}

export default function MechanismGraph({
  molecule,
  mechanismData,
  onStepSelect,
  onHighlightAtoms
}: MechanismGraphProps) {
  const [selectedStep, setSelectedStep] = useState<number | null>(null)
  const selectAtom = useMoleculeStore((s) => s.selectAtom)
  const loadMolecule = useLabStore((s) => s.loadMolecule)

  if (!mechanismData || !mechanismData.steps || mechanismData.steps.length === 0) {
    return (
      <div className="p-4 text-sm text-gray-500 text-center">
        No mechanism data available. Run mechanism prediction first.
      </div>
    )
  }

  const handleStepClick = (step: number) => {
    setSelectedStep(step)
    if (onStepSelect) {
      onStepSelect(step)
    }
  }

  const getEnergyColor = (energy: number, minEnergy: number, maxEnergy: number) => {
    if (maxEnergy === minEnergy) return 'bg-gray-100'
    const ratio = (energy - minEnergy) / (maxEnergy - minEnergy)
    if (ratio < 0.3) return 'bg-green-100 text-green-800'
    if (ratio < 0.7) return 'bg-yellow-100 text-yellow-800'
    return 'bg-red-100 text-red-800'
  }

  const energies = mechanismData.steps.map(s => s.energy)
  const minEnergy = Math.min(...energies)
  const maxEnergy = Math.max(...energies)

  return (
    <div className="space-y-4">
      <div className="flex items-center justify-between mb-4">
        <h4 className="text-sm font-semibold text-gray-900">Reaction Mechanism</h4>
        <div className="text-xs text-gray-600">
          ΔE = {mechanismData.reaction_energy.toFixed(2)} kcal/mol
        </div>
      </div>

      {/* Energy Profile */}
      <div className="p-3 bg-gray-50 rounded-lg">
        <div className="text-xs text-gray-600 mb-2">Energy Profile</div>
        <div className="relative h-32">
          {mechanismData.steps.map((step, idx) => {
            const x = (idx / (mechanismData.steps.length - 1)) * 100
            const y = 100 - ((step.energy - minEnergy) / (maxEnergy - minEnergy)) * 100
            const nextStep = mechanismData.steps[idx + 1]
            
            return (
              <React.Fragment key={step.step}>
                <div
                  className="absolute w-3 h-3 bg-blue-600 rounded-full border-2 border-white cursor-pointer hover:scale-125 transition-transform"
                  style={{ left: `${x}%`, top: `${y}%`, transform: 'translate(-50%, -50%)' }}
                  onClick={() => handleStepClick(step.step)}
                  title={`Step ${step.step}: ${step.energy.toFixed(2)} kcal/mol`}
                />
                {nextStep && (
                  <svg className="absolute inset-0 w-full h-full pointer-events-none">
                    <line
                      x1={`${x}%`}
                      y1={`${y}%`}
                      x2={`${((idx + 1) / (mechanismData.steps.length - 1)) * 100}%`}
                      y2={`${100 - ((nextStep.energy - minEnergy) / (maxEnergy - minEnergy)) * 100}%`}
                      stroke="#3b82f6"
                      strokeWidth="2"
                    />
                  </svg>
                )}
              </React.Fragment>
            )
          })}
        </div>
        <div className="flex justify-between text-xs text-gray-500 mt-2">
          <span>Reactant</span>
          <span>Product</span>
        </div>
      </div>

      {/* Steps List */}
      <div className="space-y-2 max-h-64 overflow-y-auto">
        {mechanismData.steps.map((step) => (
          <div
            key={step.step}
            className={`p-3 border rounded-lg cursor-pointer transition-all ${
              selectedStep === step.step
                ? 'border-blue-500 bg-blue-50'
                : 'border-gray-200 hover:border-gray-300 hover:bg-gray-50'
            }`}
            onClick={() => handleStepClick(step.step)}
          >
            <div className="flex items-start justify-between mb-2">
              <div className="flex-1">
                <div className="text-sm font-semibold text-gray-900">
                  Step {step.step}: {step.description}
                </div>
                <div className="text-xs text-gray-600 mt-1">
                  Energy: {step.energy.toFixed(2)} kcal/mol
                </div>
              </div>
              <div className={`px-2 py-1 text-xs rounded ${getEnergyColor(step.energy, minEnergy, maxEnergy)}`}>
                {step.energy > 0 ? '+' : ''}{step.energy.toFixed(1)}
              </div>
            </div>
            {step.intermediate && (
              <button
                onClick={(e) => {
                  e.stopPropagation()
                  // Load intermediate into viewer
                  if (step.intermediate && step.intermediate.atoms) {
                    loadMolecule(step.intermediate)
                  }
                  if (onStepSelect) {
                    onStepSelect(step.step)
                  }
                }}
                className="text-xs text-blue-600 hover:text-blue-800"
              >
                View Intermediate →
              </button>
            )}
          </div>
        ))}
      </div>

      {/* Activation Energy Info */}
      {mechanismData.activation_energy > 0 && (
        <div className="p-3 bg-yellow-50 border border-yellow-200 rounded-lg">
          <div className="text-xs font-semibold text-yellow-800 mb-1">
            Activation Energy
          </div>
          <div className="text-sm text-yellow-900">
            {mechanismData.activation_energy.toFixed(2)} kcal/mol
          </div>
        </div>
      )}
    </div>
  )
}

