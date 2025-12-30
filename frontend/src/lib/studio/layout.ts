
import type { StudioMode } from '../../types/studio';

export interface StudioLayout {
    showOptimizationPanel: boolean;
    showSimulationTimeline: boolean;
    canvasLayout: 'single' | 'split';
}

export function useStudioLayout(mode: StudioMode): StudioLayout {
    switch (mode) {
        case 'design':
            return {
                showOptimizationPanel: false,
                showSimulationTimeline: false,
                canvasLayout: 'single'
            };
        case 'optimize':
            return {
                showOptimizationPanel: true,
                showSimulationTimeline: false,
                canvasLayout: 'split'
            };
        case 'simulate':
            return {
                showOptimizationPanel: false,
                showSimulationTimeline: true,
                canvasLayout: 'single'
            };
        default:
            return {
                showOptimizationPanel: false,
                showSimulationTimeline: false,
                canvasLayout: 'single'
            };
    }
}
