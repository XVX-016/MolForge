/**
 * vitest.setup.ts
 * Global test setup: robust mocks for React/Three/React-Three-Fiber/Drei and Zustand store stubbing.
 *
 * This file:
 *  - provides lightweight DOM-safe mocks for React Three Fiber and Drei
 *  - exposes fake stores on globalThis for tests that want to opt into them
 *  - keeps setup side effects minimal so logic tests can still use real stores
 */

import { vi } from 'vitest';
import React from 'react';
import '@testing-library/jest-dom';

type PropsWithChildren = {
  children?: React.ReactNode;
};

type DivProps = PropsWithChildren & Record<string, unknown>;
type FakeStoreState = Record<string, unknown>;
type FakeStorePatch =
  | Partial<FakeStoreState>
  | ((state: FakeStoreState) => Partial<FakeStoreState>);
type MockStore = FakeStoreState & {
  getState: () => FakeStoreState;
  setState: (patch: FakeStorePatch) => void;
};

// -------------------- THREE / R3F / DREI lightweight mocks --------------------
vi.mock('@react-three/fiber', async () => {
  const React = await import('react');

  return {
    Canvas: ({ children, ...rest }: DivProps) => {
      return React.createElement('div', { 'data-testid': 'canvas', ...rest }, children);
    },
    useFrame: vi.fn((callback?: (state: { clock: { elapsedTime: number } }) => void) => {
      if (typeof callback === 'function') {
        callback({ clock: { elapsedTime: 0 } });
      }
    }),
    useThree: () => ({ camera: {}, gl: {} }),
    extend: () => {},
  };
});

vi.mock('@react-three/drei', async () => {
  const React = await import('react');

  return {
    OrbitControls: (props: Record<string, unknown>) =>
      React.createElement('div', { 'data-testid': 'orbit-controls', ...props }, null),
    Html: (props: DivProps) =>
      React.createElement('div', { 'data-testid': 'html', ...props }, null),
    Outlines: ({ color, thickness }: { color: number; thickness: number }) =>
      React.createElement('div', {
        'data-testid': 'outlines',
        'data-color': color,
        'data-thickness': thickness,
      }),
    ContactShadows: (props: Record<string, unknown>) =>
      React.createElement('div', { 'data-testid': 'contact-shadows', ...props }, null),
    Environment: (props: Record<string, unknown>) =>
      React.createElement('div', { 'data-testid': 'environment', ...props }, null),
    EffectComposer: ({ children }: PropsWithChildren) =>
      React.createElement('div', { 'data-testid': 'effect-composer' }, children),
    Bloom: (props: Record<string, unknown>) =>
      React.createElement('div', { 'data-testid': 'bloom', ...props }, null),
    ChromaticAberration: (props: Record<string, unknown>) =>
      React.createElement('div', { 'data-testid': 'chromatic-aberration', ...props }, null),
  };
});

vi.mock('three', async () => {
  class MockVector3 {
    x: number;
    y: number;
    z: number;

    constructor(x = 0, y = 0, z = 0) {
      this.x = x;
      this.y = y;
      this.z = z;
    }

    set(x: number, y: number, z: number) {
      this.x = x;
      this.y = y;
      this.z = z;
      return this;
    }

    clone() {
      return new MockVector3(this.x, this.y, this.z);
    }

    subVectors(a: MockVector3, b: MockVector3) {
      this.x = a.x - b.x;
      this.y = a.y - b.y;
      this.z = a.z - b.z;
      return this;
    }

    addVectors(a: MockVector3, b: MockVector3) {
      this.x = a.x + b.x;
      this.y = a.y + b.y;
      this.z = a.z + b.z;
      return this;
    }

    multiplyScalar(s: number) {
      this.x *= s;
      this.y *= s;
      this.z *= s;
      return this;
    }

    length() {
      return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
    }

    normalize() {
      const len = this.length();
      if (len > 0) {
        this.x /= len;
        this.y /= len;
        this.z /= len;
      }
      return this;
    }
  }

  class MockQuaternion {
    x: number;
    y: number;
    z: number;
    w: number;

    constructor(x = 0, y = 0, z = 0, w = 1) {
      this.x = x;
      this.y = y;
      this.z = z;
      this.w = w;
    }

    setFromUnitVectors() {
      return this;
    }
  }

  class MockColor {}
  class MockScene {}
  class MockGroup {}
  class MockMesh {}
  class MockPerspectiveCamera {}
  class MockWebGLRenderer {}

  return {
    Vector3: MockVector3,
    Quaternion: MockQuaternion,
    Color: MockColor,
    Scene: MockScene,
    Group: MockGroup,
    Mesh: MockMesh,
    PerspectiveCamera: MockPerspectiveCamera,
    WebGLRenderer: MockWebGLRenderer,
    MeshStandardMaterial: function MeshStandardMaterial() {
      return {};
    },
    MeshPhysicalMaterial: function MeshPhysicalMaterial() {
      return {};
    },
    SphereGeometry: function SphereGeometry() {
      return {};
    },
    CylinderGeometry: function CylinderGeometry() {
      return {};
    },
  };
});

function makeFakeStore(template: FakeStoreState): MockStore {
  const store: MockStore = { ...template } as MockStore;

  store.getState = () => store;
  store.setState = (patch: FakeStorePatch) => {
    const next = typeof patch === 'function' ? patch(store.getState()) : patch;
    Object.assign(store, next);
  };

  Object.keys(store).forEach((key) => {
    const value = store[key];
    if (typeof value === 'function' && !('_isMockFunction' in value)) {
      store[key] = vi.fn(value);
    }
  });

  return store;
}

const moleculeTemplate: FakeStoreState = {
  currentMolecule: null,
  atoms: [],
  bonds: [],
  autoBond: true,
  selectedAtomId: null,
  selectedBondId: null,
  tool: 'select',
  loadingState: 'idle',
  backendPredictions: null,
  error: null,
  atomToAdd: null,
  currentBondOrder: 1,
  addAtom: vi.fn(),
  removeAtom: vi.fn(),
  addBond: vi.fn(),
  removeBond: vi.fn(),
  selectAtom: vi.fn(),
  selectBond: vi.fn(),
  setMolecule: vi.fn((molecule: unknown) => {
    moleculeTemplate.currentMolecule = molecule;
  }),
  reset: vi.fn(() => {
    moleculeTemplate.currentMolecule = null;
    moleculeTemplate.selectedAtomId = null;
    moleculeTemplate.selectedBondId = null;
    moleculeTemplate.atoms = [];
    moleculeTemplate.bonds = [];
  }),
  setTool: vi.fn(),
  setAtomToAdd: vi.fn(),
  setBondOrder: vi.fn(),
};

const historyTemplate: FakeStoreState = {
  undoStack: {
    push: vi.fn(),
    undo: vi.fn(),
    redo: vi.fn(),
    clear: vi.fn(),
    canUndo: false,
    canRedo: false,
  },
  canUndo: false,
  canRedo: false,
};

const profileTemplate: FakeStoreState = {
  name: 'Researcher',
  avatarUrl: null,
  setName: vi.fn(),
  setAvatarUrl: vi.fn(),
  loadFromStorage: vi.fn(),
};

const fakeMoleculeStore = makeFakeStore(moleculeTemplate);
const fakeHistoryStore = makeFakeStore(historyTemplate);
const fakeProfileStore = makeFakeStore(profileTemplate);

(
  globalThis as typeof globalThis & {
    __FAKE_STORES__?: {
      moleculeStore: MockStore;
      historyStore: MockStore;
      profileStore: MockStore;
    };
  }
).__FAKE_STORES__ = {
  moleculeStore: fakeMoleculeStore,
  historyStore: fakeHistoryStore,
  profileStore: fakeProfileStore,
};
