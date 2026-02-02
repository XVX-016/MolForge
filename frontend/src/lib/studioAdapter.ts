
import type { MoleculeGraph } from '../types/molecule';

export interface RenderableAtom {
    id: string;
    element: string;
    position: [number, number, number];
}

export interface RenderableBond {
    id: string;
    from: [number, number, number];
    to: [number, number, number];
    fromAtomId: string;
    toAtomId: string;
    order: number;
}

/**
 * Converts strict MoleculeGraph to Renderable format for 3D Scene
 */
export function moleculeToRenderable(molecule: MoleculeGraph | null): {
    atoms: RenderableAtom[];
    bonds: RenderableBond[];
} {
    if (!molecule) {
        return { atoms: [], bonds: [] };
    }

    const atoms: RenderableAtom[] = molecule.atoms.map((atom) => ({
        id: atom.id,
        element: atom.element,
        position: atom.position,
    }));

    // Create a quick lookup for atom positions
    const atomMap = new Map<string, [number, number, number]>();
    molecule.atoms.forEach(a => atomMap.set(a.id, a.position));

    const bonds: RenderableBond[] = molecule.bonds.map((bond) => {
        const fromPos = atomMap.get(bond.from) || [0, 0, 0];
        const toPos = atomMap.get(bond.to) || [0, 0, 0];

        return {
            id: bond.id,
            from: fromPos as [number, number, number],
            to: toPos as [number, number, number],
            fromAtomId: bond.from,
            toAtomId: bond.to,
            order: bond.order,
        };
    }).filter(b => b.from[0] !== 0 || b.from[1] !== 0);

    return { atoms, bonds };
}

/**
 * Serializes MoleculeGraph to JSON for AI Context
 */
export function moleculeToJSON(molecule: MoleculeGraph | null): string {
    if (!molecule) return "";

    // We can simplify or strictly dump the object
    return JSON.stringify({
        atoms: molecule.atoms.map(a => ({ id: a.id, element: a.element, position: a.position })),
        bonds: molecule.bonds.map(b => ({ from: b.from, to: b.to, order: b.order })),
    });
}
