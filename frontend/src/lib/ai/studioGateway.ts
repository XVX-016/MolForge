
import { apiClient } from '../../api/api';
import type { StudioMessage, StudioMode } from '../../types/studio';

// Internal mapping of visible modes to internal AI personas
const MODE_TO_PERSONA: Record<StudioMode, string> = {
    'design': 'alina',    // Builder
    'optimize': 'ionescu', // Optimizer
    'simulate': 'solis',   // Physical Chemist
};

const PERSONA_TO_SKILL: Record<string, string> = {
    'alina': 'molecule-builder',
    'ionescu': 'property-analyst',
    'solis': 'reaction-simulator',
};

export interface ChatRequest {
    prompt: string;
    mode: StudioMode;
    moleculeContext?: string; // JSON string of molecule state
    history?: StudioMessage[];
}

export interface ChatResponse {
    response: string;
    blocked: boolean;
    error?: string;
}

/**
 * Studio AI Gateway
 * 
 * Centralizes all AI interactions for MolForge Studio.
 * Routes "Modes" to internal "Mentors" invisibly to the user.
 */
export const StudioGateway = {
    chat: async (request: ChatRequest): Promise<ChatResponse> => {
        const { prompt, mode, moleculeContext } = request;

        // 1. Resolve Mode -> Persona -> Skill
        const internalPersonaId = MODE_TO_PERSONA[mode] || 'alina';
        const mentorSkill = PERSONA_TO_SKILL[internalPersonaId];

        try {
            // 2. Wrap prompt with context if available
            let finalPrompt = prompt;
            if (moleculeContext) {
                finalPrompt = `[CONTEXT] Current Molecule JSON: ${moleculeContext} [/CONTEXT] \n\n User Request: ${prompt}`;
            }

            // 3. Call Backend
            // We use the existing mentor routing endpoint but control the persona ID internally
            const response = await apiClient.post('/api/mentor/chat', {
                prompt: finalPrompt,
                mentor_skill: mentorSkill,
                // In future, we can send 'mode' explicitly if backend supports it
            });

            return {
                response: response.data.response,
                blocked: response.data.blocked || false,
            };
        } catch (error: any) {
            console.error('StudioGateway Chat Error:', error);

            const errorMessage = error.response?.data?.detail || error.message || 'Failed to communicate with Studio AI';

            return {
                response: '',
                blocked: false,
                error: errorMessage,
            };
        }
    }
};
