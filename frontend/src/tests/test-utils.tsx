import React from 'react';
import { render } from '@testing-library/react';
import { MemoryRouter } from 'react-router-dom';

/**
 * renderWithProviders wraps the UI in MemoryRouter + Suspense
 * so tests don't need to duplicate wrappers.
 */
export function renderWithProviders(ui: React.ReactElement, opts = {}) {
  return render(
    <MemoryRouter>
      <React.Suspense fallback={null}>{ui}</React.Suspense>
    </MemoryRouter>,
    opts
  );
}


