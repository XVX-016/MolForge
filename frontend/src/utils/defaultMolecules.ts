
import type { MoleculeGraph } from '../types/molecule';

export const BENZENE: MoleculeGraph = {
    atoms: [
        { id: "a1", element: "C", position: [1.396, 0, 0], aromatic: true },
        { id: "a2", element: "C", position: [0.698, 1.209, 0], aromatic: true },
        { id: "a3", element: "C", position: [-0.698, 1.209, 0], aromatic: true },
        { id: "a4", element: "C", position: [-1.396, 0, 0], aromatic: true },
        { id: "a5", element: "C", position: [-0.698, -1.209, 0], aromatic: true },
        { id: "a6", element: "C", position: [0.698, -1.209, 0], aromatic: true },
        // Hydrogens
        { id: "h1", element: "H", position: [2.48, 0, 0] },
        { id: "h2", element: "H", position: [1.24, 2.15, 0] },
        { id: "h3", element: "H", position: [-1.24, 2.15, 0] },
        { id: "h4", element: "H", position: [-2.48, 0, 0] },
        { id: "h5", element: "H", position: [-1.24, -2.15, 0] },
        { id: "h6", element: "H", position: [1.24, -2.15, 0] },
    ],
    bonds: [
        { id: "b1", from: "a1", to: "a2", order: 2, aromatic: true },
        { id: "b2", from: "a2", to: "a3", order: 1, aromatic: true },
        { id: "b3", from: "a3", to: "a4", order: 2, aromatic: true },
        { id: "b4", from: "a4", to: "a5", order: 1, aromatic: true },
        { id: "b5", from: "a5", to: "a6", order: 2, aromatic: true },
        { id: "b6", from: "a6", to: "a1", order: 1, aromatic: true },
        // C-H bonds
        { id: "b7", from: "a1", to: "h1", order: 1 },
        { id: "b8", from: "a2", to: "h2", order: 1 },
        { id: "b9", from: "a3", to: "h3", order: 1 },
        { id: "b10", from: "a4", to: "h4", order: 1 },
        { id: "b11", from: "a5", to: "h5", order: 1 },
        { id: "b12", from: "a6", to: "h6", order: 1 },
    ],
};

export const CYCLOHEXANE: MoleculeGraph = {
    atoms: [
        { id: "c1", element: "C", position: [1.54, 0, 0] },
        { id: "c2", element: "C", position: [0.77, 1.33, 0] },
        { id: "c3", element: "C", position: [-0.77, 1.33, 0] },
        { id: "c4", element: "C", position: [-1.54, 0, 0] },
        { id: "c5", element: "C", position: [-0.77, -1.33, 0] },
        { id: "c6", element: "C", position: [0.77, -1.33, 0] },
    ],
    bonds: [
        { id: "bc1", from: "c1", to: "c2", order: 1 },
        { id: "bc2", from: "c2", to: "c3", order: 1 },
        { id: "bc3", from: "c3", to: "c4", order: 1 },
        { id: "bc4", from: "c4", to: "c5", order: 1 },
        { id: "bc5", from: "c5", to: "c6", order: 1 },
        { id: "bc6", from: "c6", to: "c1", order: 1 },
    ]
};

export const CYCLOPROPANE: MoleculeGraph = {
    atoms: [
        { id: "p1", element: "C", position: [1, -0.5, 0] },
        { id: "p2", element: "C", position: [-1, -0.5, 0] },
        { id: "p3", element: "C", position: [0, 1, 0] },
    ],
    bonds: [
        { id: "bp1", from: "p1", to: "p2", order: 1 },
        { id: "bp2", from: "p2", to: "p3", order: 1 },
        { id: "bp3", from: "p3", to: "p1", order: 1 },
    ]
};

export const CHAIN_4: MoleculeGraph = {
    atoms: [
        { id: "ch1", element: "C", position: [-2.0, 0, 0] },
        { id: "ch2", element: "C", position: [-0.6, 0.5, 0] },
        { id: "ch3", element: "C", position: [0.8, -0.5, 0] },
        { id: "ch4", element: "C", position: [2.2, 0, 0] },
    ],
    bonds: [
        { id: "bch1", from: "ch1", to: "ch2", order: 1 },
        { id: "bch2", from: "ch2", to: "ch3", order: 1 },
        { id: "bch3", from: "ch3", to: "ch4", order: 1 },
    ]
};
