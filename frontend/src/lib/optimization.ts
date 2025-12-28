
import type { MoleculeGraph } from '../types/molecule';
import { addAtom, addBond } from './mutations';

export interface OptimizationSuggestion {
    id: string;
    title: string;
    description: string;
    type: 'structure' | 'property';
    impact: {
        weightDelta: number;
        logPDelta: number;
    };
    apply: (graph: MoleculeGraph) => MoleculeGraph;
}

/**
 * Strategy: Add a Methyl Group (-CH3) to a random available position
 */
const strategyMethylate: OptimizationSuggestion = {
    id: 'strat-methylate',
    title: 'Methylate',
    description: 'Add a Methyl group (-CH3) to increase lipophilicity.',
    type: 'structure',
    impact: {
        weightDelta: 15.03, // C + 3H
        logPDelta: 0.5      // Approx
    },
    apply: (graph: MoleculeGraph) => {
        // Simple heuristic: Find a Carbon with < 4 bonds (implicit) or just pick the last Carbon
        // For this demo, we append to the last Carbon to make a chain or branch
        const carbons = graph.atoms.filter(a => a.element === 'C');
        if (carbons.length === 0) return graph;

        const target = carbons[carbons.length - 1]; // Naive target

        // Calculate new position (offset)
        const newPos: [number, number, number] = [
            target.position[0] + 1.5,
            target.position[1] + 1.0,
            target.position[2]
        ];

        // 1. Add Carbon
        let newGraph = addAtom(graph, 'C', newPos);
        const newCarbonId = newGraph.atoms[newGraph.atoms.length - 1].id;

        // 2. Bond it to target
        newGraph = addBond(newGraph, target.id, newCarbonId, 1);

        return newGraph;
    }
};

/**
 * Strategy: Add a Hydroxyl Group (-OH)
 */
const strategyHydroxylate: OptimizationSuggestion = {
    id: 'strat-hydroxylate',
    title: 'Hydroxylate',
    description: 'Add a Hydroxyl group (-OH) to increase solubility.',
    type: 'structure',
    impact: {
        weightDelta: 17.01, // O + H
        logPDelta: -0.7 // Hydrophilic
    },
    apply: (graph: MoleculeGraph) => {
        const carbons = graph.atoms.filter(a => a.element === 'C');
        if (carbons.length === 0) return graph;

        const target = carbons[0]; // Pick first Carbon for variety

        // Offset
        const newPos: [number, number, number] = [
            target.position[0] - 1.5,
            target.position[1] - 1.0,
            target.position[2]
        ];

        // 1. Add Oxygen
        let newGraph = addAtom(graph, 'O', newPos);
        const newOxygenId = newGraph.atoms[newGraph.atoms.length - 1].id;

        // 2. Bond to target
        newGraph = addBond(newGraph, target.id, newOxygenId, 1);

        return newGraph;
    }
};
export const AvailableStrategies: OptimizationSuggestion[] = [
    strategyMethylate,
    strategyHydroxylate,
];
export function generateSuggestions(graph: MoleculeGraph): OptimizationSuggestion[] {
    return AvailableStrategies;
}
