import '../../tests/setupMocks'
import { describe, it, expect } from 'vitest'
import React, { Suspense } from 'react'
import { render } from '@testing-library/react'
import { MemoryRouter } from 'react-router-dom'
import MoleculeViewer from '../MoleculeViewer'

describe('MoleculeViewer', () => {
  it('renders without crashing', () => {
    const { container } = render(
      <MemoryRouter>
        <Suspense fallback={null}>
          <MoleculeViewer />
        </Suspense>
      </MemoryRouter>
    )
    expect(container).toBeTruthy()
  })
})

