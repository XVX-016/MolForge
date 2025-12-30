
import { useState, useCallback } from 'react';
import { apiClient } from '../../api/api';
import type { StudioMessage, StudioMode } from '../../types/studio';
import { useStudioStore } from '../../store/studioStore';
import { moleculeToJSON } from '../studioAdapter';
import type { MoleculeGraph } from '../../types/molecule';

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
 * Studio AI Gateway (Static Helper)
 */
export const StudioGateway = {
    chat: async (request: ChatRequest): Promise<ChatResponse> => {
        const { prompt, mode, moleculeContext } = request;

        const internalPersonaId = MODE_TO_PERSONA[mode] || 'alina';
        const mentorSkill = PERSONA_TO_SKILL[internalPersonaId];

        try {
            let finalPrompt = prompt;
            if (moleculeContext) {
                finalPrompt = `[CONTEXT] Current Molecule JSON: ${moleculeContext} [/CONTEXT] \n\n User Request: ${prompt}`;
            }

            const response = await apiClient.post('/api/mentor/chat', {
                prompt: finalPrompt,
                mentor_skill: mentorSkill,
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

/**
 * Hook for consuming AI services in the Studio
 */
export function useStudioGateway() {
    const [isProcessing, setIsProcessing] = useState(false);
    const { addMessage } = useStudioStore();

    const sendCommand = useCallback(async (command: string, mode: StudioMode, molecule: MoleculeGraph | null) => {
        setIsProcessing(true);

        // 1. Prepare Context
        const context = molecule ? moleculeToJSON(molecule) : undefined;

        // 2. Call Gateway
        const result = await StudioGateway.chat({
            prompt: command,
            mode,
            moleculeContext: context
        });

        // 3. Handle Response
        if (result.error) {
            addMessage({
                id: Date.now().toString(),
                role: 'mentor',
                content: `Error: ${result.error}`,
                timestamp: new Date()
            });
        } else {
            addMessage({
                id: Date.now().toString(),
                role: 'mentor',
                content: result.response,
                timestamp: new Date()
            });
        }

        setIsProcessing(false);
    }, [addMessage]);

    return {
        isProcessing,
        sendCommand
    };
}
