
import type { MoleculeGraph } from '../types/molecule';

/**
 * Exports the full MoleculeGraph schema to a formatted JSON string.
 */
export function exportToJSON(graph: MoleculeGraph): string {
    return JSON.stringify(graph, null, 2);
}

/**
 * Generates a basic SMILES string for the molecule.
 * NOTE: This is a simplified generator. It handles linear chains and basic branching.
 * Complex rings and aromaticity detection are non-trivial and simplified here.
 */
export function exportToSMILES(graph: MoleculeGraph): string {
    if (graph.atoms.length === 0) return '';

    // 1. Build Adjacency Map
    const adj: Record<string, string[]> = {};
    graph.atoms.forEach(a => adj[a.id] = []);
    graph.bonds.forEach(b => {
        adj[b.from].push(b.to);
        adj[b.to].push(b.from);
    });

    const visited = new Set<string>();
    let smiles = '';

    // Depth-First Search
    function dfs(atomId: string) {
        visited.add(atomId);
        const atom = graph.atoms.find(a => a.id === atomId);
        if (!atom) return;

        smiles += atom.element;

        const neighbors = adj[atomId];
        const unvisitedNeighbors = neighbors.filter(nid => !visited.has(nid));

        // Handle Branching: All but the last neighbor get parentheses
        for (let i = 0; i < unvisitedNeighbors.length; i++) {
            const nextId = unvisitedNeighbors[i];
            const isLast = i === unvisitedNeighbors.length - 1;

            if (!isLast) smiles += '(';
            dfs(nextId);
            if (!isLast) smiles += ')';
        }
    }

    // Start from the first atom (arbitrary root)
    // A robust system would pick a better root or handle fragments.
    if (graph.atoms.length > 0) {
        dfs(graph.atoms[0].id);
    }

    return smiles;
}

/**
 * Triggers a browser download for the given content.
 */
export function downloadFile(filename: string, content: string, contentType: string = 'text/plain') {
    const blob = new Blob([content], { type: contentType });
    const url = URL.createObjectURL(blob);

    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();

    document.body.removeChild(a);
    URL.revokeObjectURL(url);
}
