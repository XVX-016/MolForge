
import type { MoleculeGraph } from '../types/molecule';

const ATOMIC_WEIGHTS: Record<string, number> = {
    H: 1.008,
    C: 12.011,
    N: 14.007,
    O: 15.999,
    F: 18.998,
    P: 30.974,
    S: 32.065,
    Cl: 35.453,
    Br: 79.904,
    I: 126.904,
};

// Simplified Crippen-like contributions for LogP estimation
// These are rough approximations for educational/demo purposes
const LOGP_CONTRIBUTIONS: Record<string, number> = {
    C: 0.2, // Generic Carbon
    H: 0.1,
    N: -0.5,
    O: -0.8, // Hydrophilic
    F: 0.1,
    Cl: 0.5,
    S: 0.3,
    Br: 0.7,
    I: 1.0
};

export function getAtomCounts(graph: MoleculeGraph): Record<string, number> {
    const counts: Record<string, number> = {};
    for (const atom of graph.atoms) {
        counts[atom.element] = (counts[atom.element] || 0) + 1;
    }
    return counts;
}

export function calculateMolecularWeight(graph: MoleculeGraph): number {
    let weight = 0;
    for (const atom of graph.atoms) {
        weight += ATOMIC_WEIGHTS[atom.element] || 0;
    }
    return Number(weight.toFixed(3));
}

export function calculateFormula(graph: MoleculeGraph): string {
    const counts = getAtomCounts(graph);

    // Hill System Order: C first, then H, then alphabetical
    let formula = '';

    if (counts['C']) {
        formula += `C${counts['C'] > 1 ? counts['C'] : ''}`;
        delete counts['C'];
    }
    if (counts['H']) {
        formula += `H${counts['H'] > 1 ? counts['H'] : ''}`;
        delete counts['H'];
    }

    const remainingElements = Object.keys(counts).sort();
    for (const el of remainingElements) {
        formula += `${el}${counts[el] > 1 ? counts[el] : ''}`;
    }

    return formula || 'Empty';
}

export function estimateLogP(graph: MoleculeGraph): number {
    let logP = 0;
    for (const atom of graph.atoms) {
        logP += LOGP_CONTRIBUTIONS[atom.element] || 0;

        // Bonus for aromatic carbons (roughly)
        if (atom.aromatic) {
            logP += 0.1;
        }
    }
    return Number(logP.toFixed(2));
}

export function calculateTPSA(graph: MoleculeGraph): number {
    // Very rough TPSA estimation based on polar atoms (N, O)
    let tpsa = 0;
    for (const atom of graph.atoms) {
        if (atom.element === 'N') tpsa += 26.02; // tertiary amine approx
        if (atom.element === 'O') tpsa += 17.07; // ether approx
    }
    return Number(tpsa.toFixed(2));
}
