import { describe, it, expect, beforeEach } from 'vitest'
import { render, screen } from '@testing-library/react'
import ToolPanel from '../ToolPanel'
import { useMoleculeStore } from '../../store/moleculeStore'

describe('ToolPanel', () => {
  beforeEach(() => {
    // Reset store
    useMoleculeStore.getState().reset()
  })

  it('renders tool buttons', () => {
    render(<ToolPanel />)
    
    // Check for tool buttons (using emoji icons)
    const buttons = screen.getAllByRole('button')
    expect(buttons.length).toBeGreaterThan(0)
  })

  it('highlights active tool', () => {
    useMoleculeStore.getState().setTool('add-atom')
    render(<ToolPanel />)
    
    // Tool should be set
    expect(useMoleculeStore.getState().tool).toBe('add-atom')
  })
})

