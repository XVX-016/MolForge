import { create } from 'zustand';
import axios from 'axios';
import type { WorkflowNode } from '../components/studio/WorkflowTimeline';

const API_BASE = 'http://localhost:8000'; // Should match backend settings

export interface Experiment {
    id: string;
    molecule_id: string;
    molecule_version_id: string;
    forcefield: string;
    created_at: string;
}

interface GeminiInsight {
    summary: string;
    key_observations: Array<{ type: string; impact: string; text: string }>;
    contradictions: Array<{ source1: string; source2: string; text: string }>;
    suggested_next_steps: Array<{ action: string; rationale: string }>;
}

interface StudioV2State {
    currentExperiment: Experiment | null;
    nodes: WorkflowNode[];
    selectedNodeId: string | null;
    compareNodeId: string | null;
    activeInsight: GeminiInsight | null;
    loading: boolean;
    polling: boolean;

    // Actions
    createExperiment: (moleculeVersionId: string) => Promise<void>;
    createExperimentFromLibrary: (libraryId: number) => Promise<void>;
    addNode: (nodeType: string, params: any, parentNodeId?: string) => Promise<void>;
    runNode: (nodeId: string) => Promise<void>;
    selectNode: (nodeId: string) => void;
    setCompareNode: (nodeId: string | null) => void;
    pollResults: () => Promise<void>;
}

export const useStudioV2Store = create<StudioV2State>((set, get) => ({
    currentExperiment: null,
    nodes: [],
    selectedNodeId: null,
    compareNodeId: null,
    activeInsight: null,
    loading: false,
    polling: false,

    createExperiment: async (moleculeVersionId: string) => {
        set({ loading: true });
        try {
            const res = await axios.post(`${API_BASE}/api/studio/v2/experiment`, null, {
                params: { molecule_version_id: moleculeVersionId }
            });
            const experiment = res.data;

            // Immediately fetch the auto-created Baseline node
            const workflowRes = await axios.get(`${API_BASE}/api/studio/v2/experiment/${experiment.id}/workflow`);

            set({
                currentExperiment: experiment,
                nodes: workflowRes.data,
                loading: false
            });
        } catch (err) {
            console.error("Failed to create experiment", err);
            set({ loading: false });
        }
    },

    createExperimentFromLibrary: async (libraryId: number) => {
        set({ loading: true });
        try {
            const res = await axios.post(`${API_BASE}/api/studio/v2/experiment/from-library/${libraryId}`);
            const experiment = res.data;

            // Immediately fetch the auto-created Baseline node
            const workflowRes = await axios.get(`${API_BASE}/api/studio/v2/experiment/${experiment.id}/workflow`);

            set({
                currentExperiment: experiment,
                nodes: workflowRes.data,
                loading: false
            });
        } catch (err) {
            console.error("Failed to create experiment from library", err);
            set({ loading: false });
        }
    },

    addNode: async (nodeType: string, params: any, parentNodeId?: string) => {
        const { currentExperiment } = get();
        if (!currentExperiment) return;

        try {
            const queryParams: Record<string, any> = { node_type: nodeType };
            if (parentNodeId) {
                queryParams.parent_node_id = parentNodeId;
            }

            const res = await axios.post(`${API_BASE}/api/studio/v2/experiment/${currentExperiment.id}/node`, params, {
                params: queryParams
            });
            set(state => ({ nodes: [...state.nodes, res.data] }));
        } catch (err) {
            console.error("Failed to add node", err);
        }
    },

    runNode: async (nodeId) => {
        try {
            await axios.post(`${API_BASE}/api/studio/v2/node/${nodeId}/run`);
            set(state => ({
                nodes: state.nodes.map(n => n.id === nodeId ? { ...n, status: 'queued' } : n)
            }));

            // Start polling if not already
            if (!get().polling) {
                get().pollResults();
            }
        } catch (err) {
            console.error("Failed to run node", err);
        }
    },

    selectNode: (nodeId) => {
        const node = get().nodes.find(n => n.id === nodeId);
        // Clear comparison if new selection is made, or keep it?
        // Let's keep it to allow changing primary while keeping comparison
        const { compareNodeId } = get();
        if (compareNodeId === nodeId) {
            // Cannot compare with self
            set({ compareNodeId: null });
        }

        // Mock insight for demo if completed
        if (node?.status === 'completed') {
            set({
                selectedNodeId: nodeId,
                activeInsight: {
                    summary: `The ${node.node_type} analysis indicates strong structural alignment but reveals potential desolvation penalties.`,
                    key_observations: [
                        { type: 'binding', impact: 'high', text: 'Affinity detected at -8.4 kcal/mol.' },
                        { type: 'stability', impact: 'medium', text: 'Binding pocket exhibits 0.88 normalized stability.' }
                    ],
                    contradictions: [],
                    suggested_next_steps: [
                        { action: 'MD', rationale: 'Validate pose stability over 50ns.' }
                    ]
                }
            });
        } else {
            set({ selectedNodeId: nodeId, activeInsight: null });
        }
    },

    setCompareNode: (nodeId) => {
        const { selectedNodeId } = get();
        if (nodeId === selectedNodeId) return; // Prevent self-comparison
        set({ compareNodeId: nodeId });
    },

    pollResults: async () => {
        const { currentExperiment } = get();
        if (!currentExperiment) return;

        set({ polling: true });

        const poll = async () => {
            try {
                const res = await axios.get(`${API_BASE}/api/studio/v2/experiment/${currentExperiment.id}/workflow`);
                const latestNodes = res.data;

                // Check if any node status changed to completed/failed
                set({ nodes: latestNodes });

                const stillRunning = latestNodes.some((n: any) => n.status === 'queued' || n.status === 'running');
                if (stillRunning) {
                    setTimeout(poll, 2000);
                } else {
                    set({ polling: false });
                }
            } catch (err) {
                console.error("Polling error", err);
                set({ polling: false });
            }
        };

        poll();
    }
}));
