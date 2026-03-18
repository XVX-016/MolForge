import type { StudioAction } from '../../types/studioActions';
import type { MoleculeGraph } from '../../types/molecule';
import type { StudioMode } from '../../types/studio';
import { apiClient } from '../../api/api';

export const STUDIO_SYSTEM_PROMPT = `
You are the MolForge AI Commander. Your role is strictly to PARSE intent from user requests into structured commands for the RDKit Chemistry Kernel.
You MUST ONLY output a JSON object. No markdown, no prose, no explanations outside the JSON.
`;

export async function processAICommand(
    input: string,
    molecule: MoleculeGraph,
    mode: StudioMode,
    analysisContext?: Record<string, unknown>
): Promise<StudioAction> {
    console.log(`[AI Control Plane] Orchestrating: "${input}" in ${mode} mode`);

    try {
        const response = await apiClient.post('/api/studio/command', {
            prompt: input,
            molecule_context: molecule,
            mode: mode,
            analysis_context: analysisContext
        });

        return response.data;
    } catch (error) {
        const errorMessage = error instanceof Error ? error.message : 'Unknown error';
        console.error('AI Control Plane Error:', error);
        return {
            type: 'NO_OP',
            reason: `Backend connection failed: ${errorMessage}`
        };
    }
}
