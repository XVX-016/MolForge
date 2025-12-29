import type { ReactionRule } from "../simulator";
import type { Bond } from "../../graph/bond";
import type { Molecule } from "../../graph/molecule";

export const Hydrogenation: ReactionRule = {
    name: "Hydrogenation",

    match(molecule: Molecule) {
        return molecule.bonds.filter(b => b.order === 2);
    },

    apply(molecule: Molecule, match: Bond) {
        const targetBond = molecule.bonds.find(b => b.id === match.id);
        if (targetBond) {
            targetBond.order = 1;
        }
        return molecule;
    }
};
