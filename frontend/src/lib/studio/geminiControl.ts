import type { StudioAction } from '../../types/studioActions';
import type { MoleculeGraph } from '../../types/molecule';
import type { StudioMode } from '../../types/studio';
import { apiClient } from '../../api/api';

export const STUDIO_SYSTEM_PROMPT = `
You are MolForge Studio AI, a rigorous molecular design orchestrator.
Use JSON only.
`;

export async function processAICommand(
    input: string,
    molecule: MoleculeGraph,
    mode: StudioMode
): Promise<StudioAction> {
    console.log(`[AI Control Plane] Orchestrating: "${input}" in ${mode} mode`);

    try {
        const response = await apiClient.post('/api/studio/command', {
            prompt: input,
            molecule_context: molecule,
            mode: mode
        });

        return response.data;
    } catch (error: any) {
        console.error('AI Control Plane Error:', error);
        return {
            type: 'NO_OP',
            reason: `Backend connection failed: ${error.message}`
        };
    }
}
