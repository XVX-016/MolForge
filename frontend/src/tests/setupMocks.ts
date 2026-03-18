import React from 'react';
import { vi } from 'vitest';

// ----- MOCK ZUSTAND STORE -----
type FakeStore = {
  currentMolecule: unknown;
  autoBond: boolean;
  selectedAtomId: string | null;
  selectedBondId: string | null;
  addAtom: ReturnType<typeof vi.fn>;
  removeAtom: ReturnType<typeof vi.fn>;
  selectAtom: ReturnType<typeof vi.fn>;
  selectBond: ReturnType<typeof vi.fn>;
  setMolecule: ReturnType<typeof vi.fn>;
  reset: ReturnType<typeof vi.fn>;
  getState: () => FakeStore;
};

type StoreSelector<T> = (store: FakeStore) => T;
type StorePatch = Partial<FakeStore> | ((store: FakeStore) => Partial<FakeStore>);

const fakeStore: FakeStore = {
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
  const useMoleculeStore = <T>(selector?: StoreSelector<T>) =>
    selector ? selector(fakeStore) : fakeStore;
  useMoleculeStore.getState = () => fakeStore;
  useMoleculeStore.setState = (partial: StorePatch) => {
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


