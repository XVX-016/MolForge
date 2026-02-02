import type { Atom, Bond } from './molecule';

export type ActionType =
    | 'CREATE_MOLECULE'
    | 'ADD_ATOM'
    | 'REMOVE_ATOM'
    | 'REPLACE_ATOM'
    | 'ADD_BOND'
    | 'REMOVE_BOND'
    | 'OPTIMIZE_GEOMETRY'
    | 'SIMULATE_REACTION'
    | 'SELECT_OPTIMIZATION_RULE'
    | 'EXPLAIN_ANALYSIS'
    | 'NO_OP';

export interface CreateMoleculeAction {
    type: 'CREATE_MOLECULE';
    payload: {
        atoms: Atom[];
        bonds: Bond[];
    };
    reason: string;
}

export interface ReplaceAtomAction {
    type: 'REPLACE_ATOM';
    payload: {
        atomId: string;
        newElement: string;
    };
    reason: string;
}

export interface AddAtomAction {
    type: 'ADD_ATOM';
    payload: {
        element: string;
        position: [number, number, number];
    };
    reason: string;
}

export interface RemoveAtomAction {
    type: 'REMOVE_ATOM';
    payload: {
        atomId: string;
    };
    reason: string;
}

export interface OptimizeGeometryAction {
    type: 'OPTIMIZE_GEOMETRY';
    reason: string;
}

export interface SimulateReactionAction {
    type: 'SIMULATE_REACTION';
    reason: string;
}

export interface AddBondAction {
    type: 'ADD_BOND';
    payload: {
        from: string;
        to: string;
        order: 1 | 2 | 3 | 'aromatic';
    };
    reason: string;
}

export interface RemoveBondAction {
    type: 'REMOVE_BOND';
    payload: {
        bondId: string;
    };
    reason: string;
}

export interface NoOpAction {
    type: 'NO_OP';
    reason: string;
}

export interface SelectOptimizationRuleAction {
    type: 'SELECT_OPTIMIZATION_RULE';
    payload: {
        rule_id: string;
    };
    reason: string;
}

export interface ExplainAnalysisAction {
    type: 'EXPLAIN_ANALYSIS';
    payload: {
        issue_code: string;
    };
    reason: string;
}

export type StudioAction =
    | CreateMoleculeAction
    | ReplaceAtomAction
    | AddAtomAction
    | RemoveAtomAction
    | AddBondAction
    | RemoveBondAction
    | OptimizeGeometryAction
    | SimulateReactionAction
    | SelectOptimizationRuleAction
    | ExplainAnalysisAction
    | NoOpAction;
