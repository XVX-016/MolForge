import type { ReactionRule } from "../simulator";

export const Esterification: ReactionRule = {
    name: "Esterification",
    match(molecule) {
        return [];
    },
    apply(molecule, match) {
        return molecule;
    }
};
