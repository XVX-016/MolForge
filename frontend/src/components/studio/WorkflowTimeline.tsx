import {
    CheckCircle2,
    Circle,
    Clock,
    AlertCircle,
    Play,
    Zap,
    Box
} from 'lucide-react';

export type WorkflowNodeState = 'created' | 'queued' | 'running' | 'completed' | 'failed' | 'invalidated';

export interface WorkflowNode {
    id: string;
    experiment_id: string;
    node_type: string;
    status: WorkflowNodeState;
    input_params: Record<string, any>;
    output_data?: Record<string, any>;
    created_at: string;
}

interface WorkflowTimelineProps {
    nodes: WorkflowNode[];
    selectedNodeId?: string;
    compareNodeId?: string;
    onNodeSelect: (nodeId: string) => void;
    onNodeCompare: (nodeId: string | null) => void;
    onRunNode: (nodeId: string) => void;
    isEmpty?: boolean;
    onAddNode: (type: string, parentId?: string) => void;
}

const GhostNode = ({ label, isBaseline = false }: { label: string, isBaseline?: boolean }) => (
    <div className={`relative flex items-start gap-3 p-3 rounded-xl border ${isBaseline ? 'bg-white border-blue-200' : 'border-dashed border-gray-200 opacity-60'
        }`}>
        <div className={`mt-0.5 rounded-full p-0.5 ${isBaseline ? 'bg-blue-100 text-blue-600' : 'bg-gray-100 text-gray-400'}`}>
            <Circle size={14} fill={isBaseline ? "currentColor" : "none"} />
        </div>
        <div>
            <span className={`text-[11px] font-black uppercase tracking-tight ${isBaseline ? 'text-black' : 'text-gray-400'}`}>
                {label}
            </span>
            {isBaseline && <p className="text-[10px] text-blue-600 font-bold mt-1">Waiting for upload...</p>}
        </div>
        {/* Connector Line */}
        <div className="absolute left-[19px] top-8 bottom-[-16px] w-[2px] bg-gray-100" />
    </div>
);

const NodeIcon = ({ type, status }: { type: string, status: WorkflowNodeState }) => {
    if (status === 'running') return <Clock className="animate-spin text-orange-500" size={18} />;
    if (status === 'completed') return <CheckCircle2 className="text-green-500" size={18} />;
    if (status === 'failed') return <AlertCircle className="text-red-500" size={18} />;

    switch (type.toLowerCase()) {
        case 'docking': return <Box className="text-blue-500" size={18} />;
        case 'md': return <Activity className="text-purple-500" size={18} />;
        case 'qsar': return <Zap className="text-yellow-500" size={18} />;
        case 'baseline': return <Circle className="text-gray-900" size={18} />;
        default: return <Circle className="text-gray-400" size={18} />;
    }
};

const Activity = ({ className, size }: { className?: string, size?: number }) => (
    <svg
        xmlns="http://www.w3.org/2000/svg"
        width={size || 24}
        height={size || 24}
        viewBox="0 0 24 24"
        fill="none"
        stroke="currentColor"
        strokeWidth="2"
        strokeLinecap="round"
        strokeLinejoin="round"
        className={className}
    >
        <path d="M22 12h-4l-3 9L9 3l-3 9H2" />
    </svg>
);

export default function WorkflowTimeline({ nodes, selectedNodeId, compareNodeId, onNodeSelect, onNodeCompare, onRunNode, onAddNode, isEmpty }: WorkflowTimelineProps) {
    if (isEmpty) {
        return (
            <div className="flex flex-col h-full bg-[#F9FAFB]">
                <div className="p-4 border-b border-[#E5E7EB] bg-white">
                    <h3 className="text-xs font-black uppercase tracking-widest text-darkGrey">Experiment Workflow</h3>
                    <p className="text-[10px] text-midGrey">Sequence of scientific analysis</p>
                </div>
                <div className="flex-1 overflow-y-auto p-4 space-y-4">
                    <GhostNode label="Baseline Molecule" isBaseline={true} />
                    <GhostNode label="Docking Screening" />
                    <GhostNode label="Molecular Dynamics" />
                    <GhostNode label="QSAR Prediction" />
                    <div className="text-center mt-6">
                        <p className="text-[10px] text-midGrey italic">Your experiment will grow here.</p>
                    </div>
                </div>
            </div>
        );
    }

    return (
        <div className="flex flex-col h-full bg-[#F9FAFB]">
            <div className="p-4 border-b border-[#E5E7EB] bg-white">
                <h3 className="text-xs font-black uppercase tracking-widest text-darkGrey">Experiment Workflow</h3>
                <p className="text-[10px] text-midGrey">Sequence of scientific analysis</p>
            </div>

            <div className="flex-1 overflow-y-auto p-4 space-y-4 relative">
                {/* Connector Line Base */}
                <div className="absolute left-[29px] top-4 bottom-4 w-[2px] bg-[#E5E7EB] z-0" />

                {nodes.map((node, index) => (
                    <div key={node.id} className="relative z-10">
                        <div
                            onClick={(e) => {
                                if (e.ctrlKey || e.metaKey) {
                                    // Toggle comparison mode
                                    if (compareNodeId === node.id) {
                                        onNodeCompare(null);
                                    } else {
                                        onNodeCompare(node.id);
                                    }
                                } else {
                                    onNodeSelect(node.id);
                                }
                            }}
                            className={`flex items-start gap-3 p-3 rounded-xl border transition-all cursor-pointer bg-white ${selectedNodeId === node.id
                                    ? 'border-blue-500 shadow-lg scale-[1.02] ring-1 ring-blue-500'
                                    : compareNodeId === node.id
                                        ? 'border-purple-500 shadow-md ring-1 ring-purple-500 bg-purple-50'
                                        : 'border-[#E5E7EB] hover:border-gray-300'
                                }`}
                        >
                            <div className="mt-0.5 bg-white rounded-full">
                                <NodeIcon type={node.node_type} status={node.status} />
                            </div>

                            <div className="flex-1 min-w-0">
                                <div className="flex justify-between items-center mb-1">
                                    <span className="text-[11px] font-black uppercase tracking-tight text-black">
                                        {node.node_type}
                                    </span>
                                    <span className={`text-[9px] font-medium px-1.5 py-0.5 rounded-full uppercase ${node.status === 'completed' ? 'bg-green-100 text-green-700' :
                                        node.status === 'running' ? 'bg-orange-100 text-orange-700' :
                                            'bg-gray-100 text-gray-700'
                                        }`}>
                                        {node.status}
                                    </span>
                                </div>
                                <p className="text-[10px] text-midGrey truncate">
                                    {Object.keys(node.input_params).length > 0 ? Object.keys(node.input_params).join(', ') : 'No params'}
                                </p>

                                {node.status === 'created' && (
                                    <button
                                        onClick={(e) => {
                                            e.stopPropagation();
                                            onRunNode(node.id);
                                        }}
                                        className="mt-2 flex items-center gap-1.5 text-[10px] font-bold text-blue-600 hover:text-blue-700"
                                    >
                                        <Play size={10} fill="currentColor" />
                                        Run Simulation
                                    </button>
                                )}
                            </div>
                        </div>
                    </div>
                ))}

                {/* Add Node Menu */}
                <div className="mt-4 pt-4 border-t border-dashed border-gray-200">
                    <p className="text-[10px] font-bold text-gray-400 mb-2 uppercase tracking-tight">Append Analysis Step</p>
                    <div className="grid grid-cols-2 gap-2">
                        <button
                            onClick={() => onAddNode('docking', selectedNodeId || undefined)}
                            className="flex flex-col items-center justify-center p-3 border border-gray-200 rounded-xl hover:bg-blue-50 hover:border-blue-200 hover:text-blue-600 transition-colors bg-white z-10"
                        >
                            <Box size={16} className="mb-1" />
                            <span className="text-[10px] font-bold uppercase">Docking</span>
                        </button>
                        <button
                            onClick={() => onAddNode('md', selectedNodeId || undefined)}
                            className="flex flex-col items-center justify-center p-3 border border-gray-200 rounded-xl hover:bg-purple-50 hover:border-purple-200 hover:text-purple-600 transition-colors bg-white z-10"
                        >
                            <Activity size={16} className="mb-1" />
                            <span className="text-[10px] font-bold uppercase">MD</span>
                        </button>
                    </div>
                </div>
            </div>
        </div>
    );
}
