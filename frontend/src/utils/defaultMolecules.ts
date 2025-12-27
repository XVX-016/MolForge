export const CH4 = {
    name: "Methane",
    atoms: [
        { id: "C1", element: "C", x: 0, y: 0, z: 0 },
        { id: "H1", element: "H", x: 1.09, y: 0, z: 0 },
        { id: "H2", element: "H", x: -0.36, y: 1.03, z: 0 },
        { id: "H3", element: "H", x: -0.36, y: -0.51, z: 0.89 },
        { id: "H4", element: "H", x: -0.36, y: -0.51, z: -0.89 },
    ],
    bonds: [
        { a: "C1", b: "H1", order: 1 },
        { a: "C1", b: "H2", order: 1 },
        { a: "C1", b: "H3", order: 1 },
        { a: "C1", b: "H4", order: 1 }
    ]
};

export const BENZENE = {
    name: "Benzene",
    atoms: [
        { id: "C1", element: "C", x: 1.396, y: 0, z: 0 },
        { id: "C2", element: "C", x: 0.698, y: 1.209, z: 0 },
        { id: "C3", element: "C", x: -0.698, y: 1.209, z: 0 },
        { id: "C4", element: "C", x: -1.396, y: 0, z: 0 },
        { id: "C5", element: "C", x: -0.698, y: -1.209, z: 0 },
        { id: "C6", element: "C", x: 0.698, y: -1.209, z: 0 },
        // Hydrogens
        { id: "H1", element: "H", x: 2.48, y: 0, z: 0 },
        { id: "H2", element: "H", x: 1.24, y: 2.15, z: 0 },
        { id: "H3", element: "H", x: -1.24, y: 2.15, z: 0 },
        { id: "H4", element: "H", x: -2.48, y: 0, z: 0 },
        { id: "H5", element: "H", x: -1.24, y: -2.15, z: 0 },
        { id: "H6", element: "H", x: 1.24, y: -2.15, z: 0 },
    ],
    bonds: [
        { a: "C1", b: "C2", order: 1 },
        { a: "C2", b: "C3", order: 2 },
        { a: "C3", b: "C4", order: 1 },
        { a: "C4", b: "C5", order: 2 },
        { a: "C5", b: "C6", order: 1 },
        { a: "C6", b: "C1", order: 2 },
        // C-H bonds
        { a: "C1", b: "H1", order: 1 },
        { a: "C2", b: "H2", order: 1 },
        { a: "C3", b: "H3", order: 1 },
        { a: "C4", b: "H4", order: 1 },
        { a: "C5", b: "H5", order: 1 },
        { a: "C6", b: "H6", order: 1 },
    ]
};

export const CYCLOHEXANE = {
    name: "Cyclohexane",
    atoms: [
        // Chair conformation approximation
        { id: "C1", element: "C", x: 1.54, y: 0, z: 0.5 },
        { id: "C2", element: "C", x: 0.77, y: 1.33, z: -0.5 },
        { id: "C3", element: "C", x: -0.77, y: 1.33, z: 0.5 },
        { id: "C4", element: "C", x: -1.54, y: 0, z: -0.5 },
        { id: "C5", element: "C", x: -0.77, y: -1.33, z: 0.5 },
        { id: "C6", element: "C", x: 0.77, y: -1.33, z: -0.5 },
    ],
    bonds: [
        { a: "C1", b: "C2", order: 1 },
        { a: "C2", b: "C3", order: 1 },
        { a: "C3", b: "C4", order: 1 },
        { a: "C4", b: "C5", order: 1 },
        { a: "C5", b: "C6", order: 1 },
        { a: "C6", b: "C1", order: 1 },
    ]
};

export const CYCLOPROPANE = {
    name: "Cyclopropane",
    atoms: [
        { id: "C1", element: "C", x: 0, y: 0.87, z: 0 },
        { id: "C2", element: "C", x: 0.75, y: -0.43, z: 0 },
        { id: "C3", element: "C", x: -0.75, y: -0.43, z: 0 },
    ],
    bonds: [
        { a: "C1", b: "C2", order: 1 },
        { a: "C2", b: "C3", order: 1 },
        { a: "C3", b: "C1", order: 1 },
    ]
};

export const CHAIN_4 = {
    name: "Butane",
    atoms: [
        { id: "C1", element: "C", x: -2.3, y: 0, z: 0 },
        { id: "C2", element: "C", x: -0.77, y: 0.5, z: 0 },
        { id: "C3", element: "C", x: 0.77, y: -0.5, z: 0 },
        { id: "C4", element: "C", x: 2.3, y: 0, z: 0 },
    ],
    bonds: [
        { a: "C1", b: "C2", order: 1 },
        { a: "C2", b: "C3", order: 1 },
        { a: "C3", b: "C4", order: 1 },
    ]
};
