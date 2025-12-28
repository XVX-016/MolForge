
import type { MoleculeGraph } from '../types/molecule';

export interface SimulationFrame {
    id: number;
    time: number; // in picoseconds (ps)
    energy: number; // in kcal/mol
    graph: MoleculeGraph;
}

export type Trajectory = SimulationFrame[];

/**
 * Generates a mock "Harmonic Vibration" trajectory.
 * Deterministic based on the input graph.
 */
export function generateVibrationTrajectory(baseGraph: MoleculeGraph, steps: number = 100): Trajectory {
    const trajectory: Trajectory = [];

    // Assign a random (but deterministic via atom index) phase and frequency to each atom
    const atomVibrations = baseGraph.atoms.map((atom, index) => ({
        phase: index * 0.5,
        freq: 1 + (index % 3) * 0.5,
        amplitude: 0.1 // Angstroms
    }));

    for (let i = 0; i < steps; i++) {
        // Clone the molecule structure
        const frameGraph: MoleculeGraph = {
            atoms: baseGraph.atoms.map((atom, idx) => {
                const vib = atomVibrations[idx];
                const t = i * 0.1; // time step

                // Perturb position
                // Simple harmonic motion around equilibrium
                const dx = Math.sin(t * vib.freq + vib.phase) * vib.amplitude;
                const dy = Math.cos(t * vib.freq + vib.phase) * vib.amplitude;
                const dz = Math.sin(t * vib.freq * 0.5) * vib.amplitude;

                return {
                    ...atom,
                    position: [
                        atom.position[0] + dx,
                        atom.position[1] + dy,
                        atom.position[2] + dz
                    ]
                };
            }),
            bonds: [...baseGraph.bonds] // Bonds don't change in simple vibration
        };

        trajectory.push({
            id: i,
            time: i * 0.1,
            energy: -100 + Math.sin(i * 0.2) * 2, // Mock energy fluctuation
            graph: frameGraph
        });
    }

    return trajectory;
}
