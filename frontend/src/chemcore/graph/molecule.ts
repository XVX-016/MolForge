import type { Atom } from "./atom";
import type { Bond } from "./bond";

export interface Molecule {
    atoms: Atom[];
    bonds: Bond[];
}
