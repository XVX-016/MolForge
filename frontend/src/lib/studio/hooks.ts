
import { useStudioStore } from '../../store/studioStore';

export interface InteractionContract {
    canEdit: boolean;
    canSelect: boolean;
    canOptimize: boolean;
    canSimulate: boolean;
    modeColor: string;
}

export function useStudioMode(): InteractionContract {
    const { mode } = useStudioStore();

    switch (mode) {
        case 'design':
            return {
                canEdit: true,
                canSelect: true,
                canOptimize: false,
                canSimulate: false,
                modeColor: '#22c55e', // Green
            };
        case 'optimize':
            return {
                canEdit: false,
                canSelect: true, // Needed for clicking suggestions/atoms
                canOptimize: true,
                canSimulate: false,
                modeColor: '#a855f7', // Purple
            };
        case 'simulate':
            return {
                canEdit: false,
                canSelect: false,
                canOptimize: false,
                canSimulate: true,
                modeColor: '#f97316', // Orange
            };
        default:
            return {
                canEdit: false,
                canSelect: false,
                canOptimize: false,
                canSimulate: false,
                modeColor: '#333333',
            };
    }
}
