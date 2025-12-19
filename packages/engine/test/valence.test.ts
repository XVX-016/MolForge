import { describe, it, expect } from 'vitest';
import { allowedAdditionalBonds, getValence } from '../src/ValenceChecker';
import { MoleculeGraph } from '../src/MoleculeGraph';
import { autoBondNewAtom } from '../src/AutoBond';

describe('ValenceChecker', () => {
  it('returns typical valences', () => {
    expect(getValence('C')).toBe(4);
    expect(getValence('N')).toBe(3);
    expect(getValence('O')).toBe(2);
    expect(getValence('H')).toBe(1);
  });

  it('computes remaining bonds', () => {
    expect(allowedAdditionalBonds('C', 0)).toBe(4);
    expect(allowedAdditionalBonds('C', 2)).toBe(2);
    expect(allowedAdditionalBonds('O', 2)).toBe(0);
    expect(allowedAdditionalBonds('H', 1)).toBe(0);
  });
});

describe('AutoBond', () => {
  it('creates bonds to nearby hydrogens around carbon', () => {
    const mol = new MoleculeGraph();
    const c = mol.addAtom({ element: 'C', position: [0, 0, 0] });
    const h1 = mol.addAtom({ element: 'H', position: [1.1, 0, 0] });
    const h2 = mol.addAtom({ element: 'H', position: [-1.1, 0, 0] });
    // Bond when new atom is carbon
    const created = autoBondNewAtom(mol, c);
    expect(created.length).toBeGreaterThan(0);
  });
});


