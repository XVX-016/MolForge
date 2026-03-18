import React from 'react'
import { renderToString } from 'react-dom/server'
import { describe, expect, it, vi } from 'vitest'

vi.mock('react-router-dom', () => ({
  Link: (props: React.ComponentProps<'a'>) => React.createElement('a', props),
  useLocation: () => ({ pathname: '/' }),
}))

import Navbar from '../Navbar'

describe('Navbar', () => {
  it('renders without crashing', () => {
    const html = renderToString(<Navbar />)
    expect(html).toContain('BioSynth AI')
  })
})
