
import type { StudioMode } from '../../types/studio';

export interface StudioLayout {
    showEditorRail: boolean;
    showActionsPanel: boolean;
    showOptimizationPanel: boolean;
    showSimulationTimeline: boolean;
    showPropertyPanel: boolean;
    canvasLayout: 'single' | 'split';
}

export function useStudioLayout(mode: StudioMode): StudioLayout {
    switch (mode) {
        case 'design':
            return {
                showEditorRail: true,
                showActionsPanel: true,
                showOptimizationPanel: false,
                showSimulationTimeline: false,
                showPropertyPanel: true,
                canvasLayout: 'single'
            };
        case 'optimize':
            return {
                showEditorRail: false,
                showActionsPanel: true,
                showOptimizationPanel: true,
                showSimulationTimeline: false,
                showPropertyPanel: true,
                canvasLayout: 'split'
            };
        case 'simulate':
            return {
                showEditorRail: false,
                showActionsPanel: false, // Playback ONLY
                showOptimizationPanel: false,
                showSimulationTimeline: true,
                showPropertyPanel: true,
                canvasLayout: 'single'
            };
        default:
            return {
                showEditorRail: true,
                showActionsPanel: true,
                showOptimizationPanel: false,
                showSimulationTimeline: false,
                showPropertyPanel: true,
                canvasLayout: 'single'
            };
    }
}
