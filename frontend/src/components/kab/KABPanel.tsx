import React, { useState } from 'react'
import { useMoleculeStore } from '../../store/moleculeStore'
import { moleculeToJSON } from '../../lib/engineAdapter'
import { MoleculeSerializer } from '@biosynth/engine'
import BindingGraph from './BindingGraph'
import AnalysisGraph from './AnalysisGraph'
import KnowledgeAlerts from './KnowledgeAlerts'

type Tab = 'binding' | 'analysis' | 'knowledge'

export default function KABPanel() {
  const molecule = useMoleculeStore((s) => s.currentMolecule)
  const [activeTab, setActiveTab] = useState<Tab>('binding')

  if (!molecule) {
    return (
      <div className="p-4 text-sm text-gray-500">
        No molecule loaded. Create or load a molecule to use KAB analysis.
      </div>
    )
  }

  return (
    <div className="flex flex-col h-full">
      {/* Tab Navigation */}
      <div className="flex border-b border-gray-200 bg-white">
        <button
          onClick={() => setActiveTab('binding')}
          className={`px-4 py-2 text-sm font-medium transition-colors ${
            activeTab === 'binding'
              ? 'text-blue-600 border-b-2 border-blue-600'
              : 'text-gray-600 hover:text-gray-900'
          }`}
        >
          Binding
        </button>
        <button
          onClick={() => setActiveTab('analysis')}
          className={`px-4 py-2 text-sm font-medium transition-colors ${
            activeTab === 'analysis'
              ? 'text-blue-600 border-b-2 border-blue-600'
              : 'text-gray-600 hover:text-gray-900'
          }`}
        >
          Analysis
        </button>
        <button
          onClick={() => setActiveTab('knowledge')}
          className={`px-4 py-2 text-sm font-medium transition-colors ${
            activeTab === 'knowledge'
              ? 'text-blue-600 border-b-2 border-blue-600'
              : 'text-gray-600 hover:text-gray-900'
          }`}
        >
          Knowledge
        </button>
      </div>

      {/* Tab Content */}
      <div className="flex-1 overflow-y-auto p-4">
        {activeTab === 'binding' && <BindingGraph molecule={molecule} />}
        {activeTab === 'analysis' && <AnalysisGraph molecule={molecule} />}
        {activeTab === 'knowledge' && <KnowledgeAlerts molecule={molecule} />}
      </div>
    </div>
  )
}

