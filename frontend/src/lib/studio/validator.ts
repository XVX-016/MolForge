import { StudioAction } from '../../types/studioActions';
import { MoleculeGraph } from '../../types/molecule';
import { StudioMode } from '../../types/studio';

export interface ValidationResult {
    valid: boolean;
    error?: string;
}

export function validateAction(
    action: StudioAction,
    molecule: MoleculeGraph,
    mode: StudioMode
): ValidationResult {
    // 1. Basic Action Presence Check
    if (!action || !action.type) {
        return { valid: false, error: 'Invalid or missing action type' };
    }

    // 2. Mode-Based Gating
    switch (mode) {
        case 'design':
            // Design mode allows most structural changes
            return validateDesignAction(action, molecule);

        case 'optimize':
            // Optimize mode ONLY allows position updates (or NO_OP)
            return validateOptimizeAction(action);

        case 'simulate':
            // Simulate mode is read-only for the structure
            return validateSimulateAction(action);

        default:
            return { valid: false, error: `Unsupported studio mode: ${mode}` };
    }
}

function validateDesignAction(action: StudioAction, molecule: MoleculeGraph): ValidationResult {
    switch (action.type) {
        case 'CREATE_MOLECULE':
            if (molecule.atoms.length > 0) {
                return { valid: false, error: 'Cannot create a new molecule while one already exists. Clear the workspace first.' };
            }
            return { valid: true };

        case 'ADD_ATOM':
        case 'REMOVE_ATOM':
        case 'REPLACE_ATOM':
        case 'ADD_BOND':
        case 'REMOVE_BOND':
        case 'NO_OP':
            return { valid: true };

        case 'OPTIMIZE_GEOMETRY':
        case 'SIMULATE_REACTION':
            return { valid: false, error: `${action.type} is not allowed in DESIGN mode. Switch modes first.` };

        default:
            return { valid: false, error: `Action ${action.type} not recognized for DESIGN mode` };
    }
}

function validateOptimizeAction(action: StudioAction): ValidationResult {
    switch (action.type) {
        case 'OPTIMIZE_GEOMETRY':
        case 'NO_OP':
            return { valid: true };

        case 'CREATE_MOLECULE':
        case 'ADD_ATOM':
        case 'REMOVE_ATOM':
        case 'REPLACE_ATOM':
        case 'ADD_BOND':
        case 'REMOVE_BOND':
        case 'SIMULATE_REACTION':
            return { valid: false, error: `Structural mutations are blocked in OPTIMIZE mode.` };

        default:
            return { valid: false, error: `Action ${action.type} not allowed in OPTIMIZE mode` };
    }
}

function validateSimulateAction(action: StudioAction): ValidationResult {
    switch (action.type) {
        case 'SIMULATE_REACTION':
        case 'NO_OP':
            return { valid: true };

        default:
            return { valid: false, error: `MolForge Studio is read-only in SIMULATE mode.` };
    }
}
