export type ElementSymbol = 'H' | 'C' | 'N' | 'O' | 'F' | 'P' | 'S' | 'Cl' | 'Br' | 'I';

const DEFAULT_VALENCE: Record<ElementSymbol, number> = {
	H: 1,
	C: 4,
	N: 3,
	O: 2,
	F: 1,
	P: 3,
	S: 2,
	Cl: 1,
	Br: 1,
	I: 1,
};

export function getValence(element: ElementSymbol): number {
	return DEFAULT_VALENCE[element] ?? 0;
}

/**
 * Given an element and the current bond order sum around it, return how many additional single bonds are allowed.
 */
export function allowedAdditionalBonds(element: ElementSymbol, currentBondOrderSum: number): number {
	const max = getValence(element);
	const remaining = Math.max(0, max - currentBondOrderSum);
	return remaining;
}


