import React from 'react'
import { renderToString } from 'react-dom/server'
import { beforeEach, describe, expect, it, vi } from 'vitest'

type ToolName = 'select' | 'add-atom'

interface MockMoleculeStore {
  tool: ToolName
  setTool: (tool: ToolName) => void
  currentBondOrder: number
  setBondOrder: (order: number) => void
  reset: () => void
}

vi.mock('react-router-dom', () => ({
  Link: (props: React.ComponentProps<'a'>) => React.createElement('a', props),
  useLocation: () => ({ pathname: '/lab' }),
}))

vi.mock('../../store/moleculeStore', () => {
  const base: MockMoleculeStore = {
    tool: 'select',
    setTool: (tool: ToolName) => {
      base.tool = tool
    },
    currentBondOrder: 1,
    setBondOrder: (order: number) => {
      base.currentBondOrder = order
    },
    reset: () => {},
  }

  const hook = <T,>(selector?: (state: MockMoleculeStore) => T) => (selector ? selector(base) : base)
  hook.getState = () => base
  hook.setState = () => {}
  return { useMoleculeStore: hook }
})

vi.mock('../../store/historyStore', () => {
  const base = { canUndo: false, canRedo: false }
  const hook = <T,>(selector?: (state: typeof base) => T) => (selector ? selector(base) : base)
  return { useHistoryStore: hook, undo: () => {}, redo: () => {} }
})

import ToolPanel from '../ToolPanel'
import { useMoleculeStore } from '../../store/moleculeStore'

describe('ToolPanel', () => {
  beforeEach(() => {
    useMoleculeStore.getState().reset()
  })

  it('renders tool panel', () => {
    const html = renderToString(<ToolPanel />)
    expect(typeof html).toBe('string')
    expect(html.length).toBeGreaterThan(0)
  })

  it('highlights active tool', () => {
    useMoleculeStore.getState().setTool('add-atom')
    renderToString(<ToolPanel />)
    expect(useMoleculeStore.getState().tool).toBe('add-atom')
  })
})
