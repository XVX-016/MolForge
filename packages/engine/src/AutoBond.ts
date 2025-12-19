import { MoleculeGraph } from './MoleculeGraph';
import type { ElementSymbol } from './ValenceChecker';
import { allowedAdditionalBonds } from './ValenceChecker';

const COVALENT_RADII: Record<ElementSymbol, number> = {
	H: 0.31,
	C: 0.76,
	N: 0.71,
	O: 0.66,
	F: 0.57,
	P: 1.07,
	S: 1.05,
	Cl: 1.02,
	Br: 1.20,
	I: 1.39,
};

const DIST_TOLERANCE = 0.45; // angstroms tolerance

function distance(a: [number, number, number], b: [number, number, number]): number {
	const dx = a[0] - b[0];
	const dy = a[1] - b[1];
	const dz = a[2] - b[2];
	return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

export function autoBondNewAtom(mol: MoleculeGraph, newAtomId: string): string[] {
	const newAtom = mol.atoms.get(newAtomId);
	if (!newAtom) return [];
	const newElement = newAtom.element as ElementSymbol;

	const bondsForAtom = (atomId: string) =>
		Array.from(mol.bonds.values()).filter((b) => b.a1 === atomId || b.a2 === atomId);

	const currentOrderSumNew = bondsForAtom(newAtomId).reduce((sum, b) => sum + b.order, 0);
	let remainingNew = allowedAdditionalBonds(newElement, currentOrderSumNew);
	if (remainingNew <= 0) return [];

	const createdBondIds: string[] = [];

	for (const [otherId, otherAtom] of mol.atoms) {
		if (otherId === newAtomId) continue;
		// skip if already bonded
		if (Array.from(mol.bonds.values()).some((b) => (b.a1 === newAtomId && b.a2 === otherId) || (b.a1 === otherId && b.a2 === newAtomId))) {
			continue;
		}
		const r1 = COVALENT_RADII[newElement] ?? 0.8;
		const r2 = COVALENT_RADII[otherAtom.element as ElementSymbol] ?? 0.8;
		const maxDist = r1 + r2 + DIST_TOLERANCE;
		const d = distance(newAtom.position, otherAtom.position);
		if (d > maxDist) continue;

		// valence check for neighbor
		const neighborOrderSum = bondsForAtom(otherId).reduce((sum, b) => sum + b.order, 0);
		const remainingNeighbor = allowedAdditionalBonds(otherAtom.element as ElementSymbol, neighborOrderSum);
		if (remainingNeighbor <= 0) continue;

		// create single bond
		const id = mol.addBond(newAtomId, otherId, 1);
		if (id) {
			createdBondIds.push(id);
			remainingNew -= 1;
			if (remainingNew <= 0) break;
		}
	}
	return createdBondIds;
}


