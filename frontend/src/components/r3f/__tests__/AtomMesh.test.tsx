import { describe, it, expect, vi } from 'vitest';
import React from 'react';
import { renderToString } from 'react-dom/server';
vi.mock('react-router-dom', () => {
  const React = require('react');
  return {
    Link: (props: any) => React.createElement('a', props),
    useLocation: () => ({ pathname: '/' }),
  };
});
vi.mock('../../store/moleculeStore', () => {
  const base = { selectedAtomId: null, tool: 'select' };
  const hook: any = (selector?: any) => (selector ? selector(base) : base);
  return { useMoleculeStore: hook };
});
import AtomMesh from '../AtomMesh';

describe('AtomMesh', () => {
  it('renders without throwing', () => {
    const html = renderToString(<AtomMesh id="a1" position={[0,0,0]} element="C" />);
    expect(typeof html).toBe('string');
  });
});


