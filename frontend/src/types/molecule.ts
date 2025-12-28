
export type AtomId = string;
export type BondId = string;

export interface Atom {
  id: AtomId;
  element: "C" | "H" | "O" | "N" | "S" | "Cl" | "F" | string;
  position: [number, number, number];
  charge?: number;
  aromatic?: boolean;
}

export interface Bond {
  id: BondId;
  from: AtomId;
  to: AtomId;
  order: 1 | 2 | 3;
  aromatic?: boolean;
}

export interface MoleculeGraph {
  atoms: Atom[];
  bonds: Bond[];
}
