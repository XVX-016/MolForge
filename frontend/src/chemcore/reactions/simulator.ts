import type { Molecule } from "../graph/molecule";

export type ReactionMatch = Record<string, unknown>;

export interface ReactionRule {
    name: string;
    // Match returns array of 'match objects' (could be Atoms, Bonds, or custom structure)
    match(molecule: Molecule): ReactionMatch[];
    apply(molecule: Molecule, match: ReactionMatch): Molecule;
}

export function simulateReaction(
    molecule: Molecule,
    rule: ReactionRule
): Molecule | null {
    const matches = rule.match(molecule);
    if (!matches.length) return null;

    // Deep clone to ensure immutability
    const cloned = JSON.parse(JSON.stringify(molecule)) as Molecule;

    // Apply to the first match for now (deterministic single step)
    return rule.apply(cloned, matches[0]);
}
