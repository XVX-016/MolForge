import type { Molecule as CoreMolecule } from "../chemcore/graph/molecule";
import type { Atom } from "../chemcore/graph/atom";
import type { Bond } from "../chemcore/graph/bond";

export type { Atom, Bond };

export interface Molecule extends CoreMolecule {
  id: string;
  name?: string;
  formula?: string;
  previewImage?: string;
  createdAt?: string;
  updatedAt?: string;
  metadata?: Record<string, any>;
  isValid?: boolean;
  qualityScore?: number;
}

export type ToolName =
  | 'select'
  | 'add-atom'
  | 'add-bond'
  | 'element'
  | 'erase'
  | 'rotate';

export type Vec3 = [number, number, number];
