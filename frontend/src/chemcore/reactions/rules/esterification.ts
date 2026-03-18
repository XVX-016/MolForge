import type { ReactionRule } from "../simulator";

export const Esterification: ReactionRule = {
    name: "Esterification",
    match(molecule) {
        void molecule;
        return [];
    },
    apply(molecule, match) {
        void match;
        return molecule;
    }
};
