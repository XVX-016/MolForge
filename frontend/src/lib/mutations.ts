
import type { MoleculeGraph, Atom, Bond, AtomId, BondId } from '../types/molecule';

/**
 * Pure function to add an atom to the graph.
 * Returns a NEW MoleculeGraph.
 */
export function addAtom(
    graph: MoleculeGraph,
    element: string,
    position: [number, number, number]
): MoleculeGraph {
    const newId: AtomId = `a-${Date.now()}-${Math.random().toString(36).substr(2, 5)}`;

    const newAtom: Atom = {
        id: newId,
        element,
        position,
    };

    return {
        ...graph,
        atoms: [...graph.atoms, newAtom],
        bonds: [...graph.bonds]
    };
}

/**
 * Pure function to delete an atom and its connected bonds.
 */
export function deleteAtom(graph: MoleculeGraph, atomId: AtomId): MoleculeGraph {
    // 1. Filter out the atom
    const newAtoms = graph.atoms.filter(a => a.id !== atomId);

    // 2. Filter out any bonds connected to this atom
    const newBonds = graph.bonds.filter(b => b.from !== atomId && b.to !== atomId);

    return {
        atoms: newAtoms,
        bonds: newBonds
    };
}

/**
 * Pure function to add a bond between two atoms.
 */
export function addBond(
    graph: MoleculeGraph,
    fromId: AtomId,
    toId: AtomId,
    order: 1 | 2 | 3
): MoleculeGraph {
    // Prevent duplicate bonds
    const exists = graph.bonds.some(b =>
        (b.from === fromId && b.to === toId) ||
        (b.from === toId && b.to === fromId)
    );

    if (exists) return graph;

    const newBond: Bond = {
        id: `b-${Date.now()}-${Math.random().toString(36).substr(2, 5)}`,
        from: fromId,
        to: toId,
        order
    };

    return {
        ...graph,
        atoms: [...graph.atoms],
        bonds: [...graph.bonds, newBond]
    };
}

/**
 * Pure function to delete a bond.
 */
export function deleteBond(graph: MoleculeGraph, bondId: BondId): MoleculeGraph {
    return {
        ...graph,
        atoms: [...graph.atoms],
        bonds: graph.bonds.filter(b => b.id !== bondId)
    };
}

/**
 * Pure function to update an atom's element.
 */
export function updateAtomElement(graph: MoleculeGraph, atomId: AtomId, newElement: string): MoleculeGraph {
    return {
        ...graph,
        atoms: graph.atoms.map(a => a.id === atomId ? { ...a, element: newElement } : a),
        bonds: [...graph.bonds]
    };
}
