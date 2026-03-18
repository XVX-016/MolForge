import axios from 'axios';
import type { Atom, Bond } from '../types/molecule';

interface MoleculeData {
    atoms: Atom[];
    bonds: Bond[];
}

export async function optimizeGeometry(molecule: MoleculeData) {
    try {
        const { data } = await axios.post("/api/ml/optimize", molecule);
        return data;
    } catch (error) {
        console.error("Geometry optimization failed:", error);
        throw error;
    }
}

export async function validateBond(molecule: MoleculeData, bond: Bond) {
    try {
        const { data } = await axios.post("/api/ml/validate-bond", { molecule, bond });
        return data;
    } catch (error) {
        console.error("Bond validation failed:", error);
        throw error;
    }
}
