import * as OCL from "openchemlib";
import { nanoid } from "nanoid";
import type { Molecule } from "../graph/molecule";
import type { Atom } from "../graph/atom";
import type { Bond } from "../graph/bond";

export function fromSmiles(smiles: string): Molecule {
    const mol = OCL.Molecule.fromSmiles(smiles);

    const atomCount = mol.getAllAtoms();
    const bondCount = mol.getAllBonds();

    const atoms: Atom[] = [];
    const bonds: Bond[] = [];
    const atomIdMap = new Map<number, string>(); // Index -> ID

    // 1. Parse Atoms
    for (let i = 0; i < atomCount; i++) {
        const id = nanoid();
        atomIdMap.set(i, id);

        atoms.push({
            id,
            element: mol.getAtomLabel(i) || 'C', // Fallback
            position: {
                x: mol.getAtomX(i),
                y: mol.getAtomY(i),
                z: mol.getAtomZ(i)
            }
        });
    }

    // 2. Parse Bonds
    for (let i = 0; i < bondCount; i++) {
        const fromIndex = mol.getBondAtom(0, i);
        const toIndex = mol.getBondAtom(1, i);
        const order = mol.getBondOrder(i); // OCL: 1, 2, 3

        bonds.push({
            id: nanoid(),
            from: atomIdMap.get(fromIndex)!,
            to: atomIdMap.get(toIndex)!,
            order: (order >= 1 && order <= 3) ? order as 1 | 2 | 3 : 1
        });
    }

    return { atoms, bonds };
}

export function toSmiles(molecule: Molecule): string {
    const mol = new OCL.Molecule(molecule.atoms.length, molecule.bonds.length);
    const atomIndexMap = new Map<string, number>();

    // Add Atoms
    molecule.atoms.forEach((atom) => {
        const i = mol.addAtom(OCL.Molecule.getAtomicNoFromLabel(atom.element));
        atomIndexMap.set(atom.id, i);
        if (atom.position) {
            mol.setAtomX(i, atom.position.x);
            mol.setAtomY(i, atom.position.y);
            mol.setAtomZ(i, atom.position.z);
        }
    });

    // Add Bonds
    molecule.bonds.forEach(bond => {
        const from = atomIndexMap.get(bond.from);
        const to = atomIndexMap.get(bond.to);
        if (from !== undefined && to !== undefined) {
            const bondIndex = mol.addBond(from, to);
            mol.setBondOrder(bondIndex, bond.order);
        }
    });

    return mol.toSmiles();
}
