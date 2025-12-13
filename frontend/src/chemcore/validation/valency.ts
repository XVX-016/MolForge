import type { Molecule } from "../graph/molecule";

export const VALENCY: Record<string, number> = {
    H: 1,
    C: 4,
    N: 3,
    O: 2,
    S: 6,
    P: 5,
    F: 1,
    Cl: 1,
    Br: 1,
    I: 1
};

export function getAtomBondOrderSum(
    molecule: Molecule,
    atomId: string
) {
    return molecule.bonds.reduce((sum, bond) => {
        if (bond.from === atomId || bond.to === atomId) {
            return sum + bond.order;
        }
        return sum;
    }, 0);
}

export function canAddBond(
    molecule: Molecule,
    atomId: string,
    newBondOrder: number
) {
    const atom = molecule.atoms.find(a => a.id === atomId);
    if (!atom) return false;

    const maxValence = VALENCY[atom.element];
    if (maxValence === undefined) return true; // Unknown element? Allow for now (or strict block?) -> Allow.

    const used = getAtomBondOrderSum(molecule, atomId);
    return used + newBondOrder <= maxValence;
}

export interface ValidationResult {
    valid: boolean;
    atomId?: string;
    max?: number;
    used?: number;
    reason?: string;
}

export function validateValency(molecule: Molecule): ValidationResult {
    const usage: Record<string, number> = {};

    for (const bond of molecule.bonds) {
        usage[bond.from] = (usage[bond.from] || 0) + bond.order;
        usage[bond.to] = (usage[bond.to] || 0) + bond.order;
    }

    for (const atom of molecule.atoms) {
        const used = usage[atom.id] || 0;
        const max = VALENCY[atom.element];

        if (max !== undefined && used > max) {
            return {
                valid: false,
                atomId: atom.id,
                max,
                used,
                reason: `Valency exceeded for ${atom.element} (Max ${max}, Used ${used})`
            };
        }
    }

    return { valid: true };
}
