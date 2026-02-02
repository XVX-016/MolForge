import Studio3DScene from './Studio3DScene';
import type { WorkflowNode } from './WorkflowTimeline';
import { useStudioV2Store } from '../../store/useStudioV2Store';
import axios from 'axios';

interface MultimodalWorkspaceProps {
    activeNode?: WorkflowNode;
    compareNode?: WorkflowNode;
    onCloseComparison?: () => void;
    _moleculeData: any;
    isEmpty?: boolean;
}

const EmptyState = () => {
    const { createExperimentFromLibrary } = useStudioV2Store();

    const handleLoadMolecule = async () => {
        try {
            // Quick hack: Fetch first molecule from library to start
            const res = await axios.get('http://localhost:8000/molecules/list?limit=1');
            if (res.data && res.data.length > 0) {
                await createExperimentFromLibrary(res.data[0].id);
            } else {
                alert("Library is empty. Please upload a molecule to the Library first.");
            }
        } catch (err) {
            console.error("Failed to load molecule", err);
            alert("Failed to connect to library service.");
        }
    };

    return (
        <div className="flex-1 flex flex-col items-center justify-center p-8 text-center relative overflow-hidden">
            {/* Ghost Grid Background */}
            <div className="absolute inset-0 z-0 opacity-[0.03]"
                style={{ backgroundImage: 'radial-gradient(circle, #000 1px, transparent 1px)', backgroundSize: '20px 20px' }}
            />

            <div className="z-10 bg-white p-8 rounded-2xl border border-[#E5E7EB] shadow-xl max-w-md w-full">
                <div className="w-16 h-16 bg-gray-100 rounded-full flex items-center justify-center mx-auto mb-6">
                    <div className="w-8 h-8 border-2 border-gray-300 rounded-lg transform rotate-45" />
                </div>

                <h3 className="text-xl font-black text-black mb-2 tracking-tight">Untitled Experiment</h3>
                <p className="text-sm text-midGrey mb-8 leading-relaxed">
                    Load a molecule to initialize the scientific kernel and begin the discovery workflow.
                </p>

                <button
                    onClick={handleLoadMolecule}
                    className="w-full py-3 bg-blue-600 hover:bg-blue-700 text-white font-bold rounded-xl transition-all shadow-lg shadow-blue-500/20 active:scale-95 mb-4"
                >
                    Load Molecule
                </button>

                <div className="flex justify-center gap-4 text-[10px] font-bold text-gray-400 uppercase tracking-widest">
                    <span className="hover:text-blue-500 cursor-pointer transition-colors">Import SDF</span>
                    <span className="w-px h-3 bg-gray-200" />
                    <span className="hover:text-blue-500 cursor-pointer transition-colors">From Library</span>
                </div>
            </div>

            {/* Disabled Metadata Table Placeholder */}
            <div className="absolute bottom-0 inset-x-0 h-12 bg-white border-t border-[#E5E7EB] flex items-center justify-center grayscale opacity-50 pointer-events-none">
                <div className="flex gap-12 text-[10px] font-mono text-gray-400">
                    <span>MW: —</span>
                    <span>LogP: —</span>
                    <span>TPSA: —</span>
                </div>
            </div>
        </div>
    );
};

const DockingView = (_props: { node: WorkflowNode }) => (
    <div className="absolute inset-0 flex flex-col">
        <div className="flex-1 relative">
            <Studio3DScene mode="optimize" molecule={null} />
        </div>
        <div className="h-[200px] border-t bg-white p-4 overflow-y-auto">
            <h4 className="text-[10px] font-black uppercase text-darkGrey mb-2">Binding Poses</h4>
            <table className="w-full text-[11px]">
                <thead>
                    <tr className="border-b text-midGrey">
                        <th className="text-left font-bold py-2">Pose ID</th>
                        <th className="text-right font-bold py-2">Affinity (kcal/mol)</th>
                        <th className="text-right font-bold py-2">RMSD</th>
                    </tr>
                </thead>
                <tbody>
                    <tr className="border-b">
                        <td className="py-2">Pose 0 (Core)</td>
                        <td className="text-right font-mono font-bold text-blue-600">-8.4</td>
                        <td className="text-right font-mono">0.00</td>
                    </tr>
                </tbody>
            </table>
        </div>
    </div>
);

const MDView = (_props: { node: WorkflowNode }) => (
    <div className="absolute inset-0 bg-[#F9FAFB] p-6 flex flex-col gap-6 overflow-y-auto">
        <h3 className="text-sm font-black text-black uppercase">Molecular Dynamics Stability Profile</h3>
        <div className="h-[300px] bg-white border rounded-2xl p-4 shadow-sm flex items-center justify-center">
            <span className="text-midGrey italic">[RMSD / RMSF Plot Component Integration]</span>
        </div>
        <div className="grid grid-cols-2 gap-4">
            <div className="bg-white p-4 border rounded-2xl shadow-sm">
                <span className="text-[10px] uppercase font-black text-midGrey">Stability Score</span>
                <p className="text-2xl font-mono font-bold text-green-600">0.88</p>
            </div>
            <div className="bg-white p-4 border rounded-2xl shadow-sm">
                <span className="text-[10px] uppercase font-black text-midGrey">Convergence</span>
                <p className="text-2xl font-mono font-bold text-blue-600">PASSED</p>
            </div>
        </div>
    </div>
);

// Helper to render specific node views
const renderNodeView = (node: WorkflowNode) => {
    switch (node.node_type) {
        case 'docking':
            return <DockingView node={node} />;
        case 'md':
            return <MDView node={node} />;
        default:
            return (
                <div className="flex-1 relative">
                    <Studio3DScene mode="optimize" molecule={null} />
                </div>
            );
    }
};

export default function MultimodalWorkspace({ activeNode, compareNode, onCloseComparison, _moleculeData, isEmpty }: MultimodalWorkspaceProps) {
    if (isEmpty) {
        return <EmptyState />;
    }

    // Silence unused variable warning until property is implemented
    void _moleculeData;

    if (!activeNode) {
        return (
            <div className="flex-1 relative">
                <Studio3DScene mode="optimize" molecule={null} />
            </div>
        );
    }

    // Comparison Mode: Split View
    if (compareNode) {
        return (
            <div className="flex-1 flex flex-row h-full relative">
                {/* Delta Badge Overlay */}
                {activeNode.node_type === compareNode.node_type && (
                    <div className="absolute top-1/2 left-1/2 -translate-x-1/2 -translate-y-1/2 z-20 bg-white border border-gray-200 shadow-lg rounded-full px-4 py-2 flex items-center gap-3">
                        <span className="text-[10px] uppercase font-black text-midGrey">Delta</span>
                        <div className="flex items-center gap-1 font-mono text-xs font-bold text-black">
                            <span>-0.2</span>
                            <span className="text-green-500 text-[9px]">(Improved)</span>
                        </div>
                    </div>
                )}

                {/* Primary Panel (Left) */}
                <div className="flex-1 flex flex-col border-r border-[#E5E7EB] relative">
                    <div className="absolute top-4 left-4 z-10 bg-blue-600 text-white text-[10px] font-bold px-2 py-1 rounded shadow-md uppercase tracking-widest">
                        Primary: {activeNode.node_type}
                    </div>
                    {renderNodeView(activeNode)}
                </div>

                {/* Comparison Panel (Right) */}
                <div className="flex-1 flex flex-col relative bg-gray-50">
                    <div className="absolute top-4 left-4 z-10 bg-purple-600 text-white text-[10px] font-bold px-2 py-1 rounded shadow-md uppercase tracking-widest flex items-center gap-2">
                        <span>Comparison: {compareNode.node_type}</span>
                        {onCloseComparison && (
                            <button
                                onClick={onCloseComparison}
                                className="bg-purple-700 hover:bg-purple-800 rounded px-1.5 py-0.5 transition-colors"
                            >
                                ✕
                            </button>
                        )}
                    </div>
                    {renderNodeView(compareNode)}
                </div>
            </div>
        );
    }

    // Single View
    return (
        <div className="flex-1 flex flex-col h-full relative">
            {renderNodeView(activeNode)}
        </div>
    );
}
