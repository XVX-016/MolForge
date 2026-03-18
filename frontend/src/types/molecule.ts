export type AtomId = string
export type BondId = string
export type ToolName = 'select' | 'add_atom' | 'bond' | 'delete' | 'move' | 'inspect'
export type Vec3 = [number, number, number]

export interface Position3D {
  x: number
  y: number
  z: number
}

export interface Atom {
  id: AtomId
  element: 'C' | 'H' | 'O' | 'N' | 'S' | 'Cl' | 'F' | string
  position: Position3D
  charge?: number
  aromatic?: boolean
}

export interface Bond {
  id: BondId
  from: AtomId
  to: AtomId
  order: 1 | 1.5 | 2 | 3
  aromatic?: boolean
}

export interface MoleculeMetadata {
  name?: string
  smiles?: string
  formula?: string
  molfile?: string
  source?: 'user' | 'public' | string
  [key: string]: unknown
}

export interface Molecule {
  id: string
  atoms: Atom[]
  bonds: Bond[]
  metadata?: MoleculeMetadata
}

export type MoleculeGraph = Molecule
