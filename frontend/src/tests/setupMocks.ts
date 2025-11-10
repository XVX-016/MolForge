import React from 'react';
import { vi } from 'vitest';

// ----- MOCK ZUSTAND STORE -----
const fakeStore: any = {
  currentMolecule: null,
  autoBond: true,
  selectedAtomId: null,
  selectedBondId: null,
  addAtom: vi.fn(),
  removeAtom: vi.fn(),
  selectAtom: vi.fn(),
  selectBond: vi.fn(),
  setMolecule: vi.fn((molecule: unknown) => {
    fakeStore.currentMolecule = molecule;
  }),
  reset: vi.fn(() => {
    fakeStore.currentMolecule = null;
    fakeStore.selectedAtomId = null;
    fakeStore.selectedBondId = null;
  }),
  getState: () => fakeStore,
};

vi.mock('../store/moleculeStore', () => {
  const useMoleculeStore: any = (selector?: (store: typeof fakeStore) => unknown) =>
    selector ? selector(fakeStore) : fakeStore;
  useMoleculeStore.getState = () => fakeStore;
  useMoleculeStore.setState = (partial: any) => {
    const next = typeof partial === 'function' ? partial(fakeStore) : partial;
    Object.assign(fakeStore, next);
  };
  return { useMoleculeStore };
});

// ----- MOCK R3F / THREE PRIMITIVES -----
vi.mock('@react-three/fiber', () => ({
  Canvas: ({ children }: { children?: React.ReactNode }) => React.createElement('div', null, children),
  useFrame: () => {},
  useThree: () => ({ camera: {}, gl: {} }),
}));

vi.mock('@react-three/drei', () => ({
  OrbitControls: () => React.createElement('div'),
  Html: ({ children }: { children?: React.ReactNode }) => React.createElement('div', null, children),
  Outlines: () => React.createElement('div'),
  ContactShadows: () => React.createElement('div'),
  Environment: () => React.createElement('div'),
}));

vi.mock('../test/r3f-mock');


