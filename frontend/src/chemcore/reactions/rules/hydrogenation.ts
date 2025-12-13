import type { ReactionRule } from "../simulator";
import type { Bond } from "../../graph/bond";
import type { Molecule } from "../../graph/molecule";

export const Hydrogenation: ReactionRule = {
    name: "Hydrogenation",

    match(molecule: Molecule) {
        // Find all double bonds
        return molecule.bonds.filter(b => b.order === 2);
    },

    apply(molecule: Molecule, match: Bond) {
        // match is the bond object from the CLONED molecule (passed by simulator? No, match found in original.)
        // Wait, simulator clones AFTER match.
        // So 'match' refers to object in ORIGINAL molecule.
        // We need to find the corresponding bond in the CLONED molecule.
        // IDs are stable.

        // Safety check if match is indeed a bond with ID
        const targetBond = molecule.bonds.find(b => b.id === match.id);
        if (targetBond) {
            targetBond.order = 1;
        }
        return molecule;
    }
};
