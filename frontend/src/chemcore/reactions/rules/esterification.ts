import type { ReactionRule } from "../simulator";

export const Esterification: ReactionRule = {
    name: "Esterification",
    match(molecule) {
        // TODO: Implement substructure matching for Carboxylic Group (R-COOH) and Alcohol (R-OH)
        // Requires Graph Substructure Search algorithim (O(N^2) or Ullmann)
        return [];
    },
    apply(molecule, match) {
        return molecule;
    }
};
