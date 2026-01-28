
import { create } from 'zustand';
import { apiClient } from '../api/api';

export type StudioStatus =
    | "IDLE"
    | "LOADING"
    | "READY"
    | "OPTIMIZING"
    | "PROPOSED"
    | "COMMITTED";

export type SelectionType = 'atom' | 'bond' | null;

export interface StudioSelection {
    type: SelectionType;
    id: string | null;
}

export interface DashboardPayload {
    baseline: {
        version_id: string;
        smiles: string;
        properties: Record<string, any>;
    };
    proposal: {
        version_id: string;
        smiles: string;
        properties: Record<string, any>;
    } | null;
    diff: {
        atoms: { added: number[]; removed: number[]; modified: number[] };
        bonds: { added: number[][]; removed: number[][] };
    };
    alerts: Array<{ id: string; severity: string; description: string }>;
    property_delta: Record<string, number>;
    radar: {
        baseline: Record<string, number>;
        proposal: Record<string, number>;
    };
    optimization_context: {
        available_rules: Array<{
            id: string;
            title: string;
            description: string;
            rationale: string;
            impact: Record<string, number>;
        }>;
    };
}

interface StudioState {
    status: StudioStatus;
    baselineVersionId: string | null;
    proposalVersionId: string | null;
    dashboard: DashboardPayload | null;
    error: string | null;
    selection: StudioSelection;

    // Actions
    loadDashboard: (baselineId: string, proposalId?: string | null) => Promise<void>;
    runAICommand: (intent: string) => Promise<void>;
    applyRule: (ruleId: string) => Promise<void>;
    acceptProposal: () => Promise<void>;
    rejectProposal: () => Promise<void>;
    setSelection: (type: SelectionType, id: string | null) => void;
}

export const useStudioStore = create<StudioState>((set, get) => ({
    status: "IDLE",
    baselineVersionId: null,
    proposalVersionId: null,
    dashboard: null,
    error: null,
    selection: { type: null, id: null },

    setSelection: (type, id) => set({ selection: { type, id } }),

    loadDashboard: async (baselineId, proposalId = null) => {
        set({ status: "LOADING", error: null });
        try {
            const response = await apiClient.post("/api/molecule/dashboard", {
                baseline_version_id: baselineId,
                proposal_version_id: proposalId
            });
            set({
                dashboard: response.data,
                baselineVersionId: baselineId,
                proposalVersionId: proposalId,
                status: proposalId ? "PROPOSED" : "READY"
            });
        } catch (err: any) {
            set({ status: "IDLE", error: err.message || "Failed to load dashboard" });
        }
    },

    runAICommand: async (intent) => {
        const { status, baselineVersionId, dashboard } = get();
        if (status !== "READY" || !baselineVersionId) return;

        set({ status: "OPTIMIZING", error: null });
        try {
            // 1. Get AI proposal
            const aiRes = await apiClient.post("/api/studio/command", {
                prompt: intent,
                molecule_context: { smiles: dashboard?.baseline.smiles },
                mode: "optimize",
                analysis_context: dashboard?.optimization_context
            });

            const action = aiRes.data;
            if (action.type === "SELECT_OPTIMIZATION_RULE") {
                await get().applyRule(action.payload.rule_id);
            } else {
                set({ status: "READY", error: action.payload?.reason || "AI found no appropriate optimization." });
            }
        } catch (err: any) {
            set({ status: "READY", error: err.message || "AI Command failed" });
        }
    },

    applyRule: async (_ruleId) => {
        const { status, baselineVersionId } = get();
        if ((status !== "READY" && status !== "OPTIMIZING") || !baselineVersionId) return;

        set({ status: "OPTIMIZING", error: null });
        try {
            // In this authoritative model, 'applying' a rule means asking the backend
            // to create an ephemeral/optimized version and then loading the dashboard with it.

            await apiClient.get(`/api/molecule/dashboard`); // Placeholder for Rule execution

            set({ status: "PROPOSED" }); // Mock transition
        } catch (err: any) {
            set({ status: "READY", error: err.message || "Rule application failed" });
        }
    },

    acceptProposal: async () => {
        const { status, proposalVersionId } = get();
        if (status !== "PROPOSED" || !proposalVersionId) return;

        set({ status: "LOADING" });
        try {
            // Transition ephemeral OPTIMIZED state to VERSIONED
            await apiClient.post(`/api/molecule/${proposalVersionId}/accept`);
            set({ status: "COMMITTED" });
            setTimeout(() => set({ status: "READY", proposalVersionId: null }), 2000);
        } catch (err: any) {
            set({ status: "PROPOSED", error: err.message || "Failed to accept proposal" });
        }
    },

    rejectProposal: async () => {
        set({ status: "READY", proposalVersionId: null, error: null });
    }
}));
