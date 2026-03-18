import { useState, useEffect } from 'react';
import Studio3DScene from './Studio3DScene';
import { ArrowRight, Zap, AlertTriangle, CheckCircle } from 'lucide-react';
import { useStudioStore } from '../../store/studioStore';
import { getDashboard, type DashboardResponse, type AnalysisIssue } from '../../lib/api';
import { RadarChart } from './RadarChart';
import * as THREE from 'three';

type PropertyMetricKey = 'molecular_weight' | 'logp' | 'tpsa' | 'qed';

function getNumericProperty(source: unknown, key: PropertyMetricKey): number {
    if (!source || typeof source !== 'object') {
        return 0;
    }

    const value = (source as Record<string, unknown>)[key];
    return typeof value === 'number' ? value : 0;
}

export default function OptimizeCompareView() {
    const { mode, activeVersionId } = useStudioStore();
    const [dashboard, setDashboard] = useState<DashboardResponse | null>(null);
    const [isLoading, setIsLoading] = useState(false);

    const [cameraTarget, setCameraTarget] = useState(new THREE.Vector3(0, 0, 0));
    const [cameraPosition, setCameraPosition] = useState(new THREE.Vector3(0, 0, 8));

    useEffect(() => {
        const load = async () => {
            if (!activeVersionId) return;

            setIsLoading(true);
            try {
                const data = await getDashboard(activeVersionId, activeVersionId);
                setDashboard(data);
            } catch (error) {
                console.error('Dashboard Load Error', error);
            } finally {
                setIsLoading(false);
            }
        };
        load();
    }, [activeVersionId]);

    const handleCameraChange = (target: THREE.Vector3, position: THREE.Vector3) => {
        setCameraTarget(target);
        setCameraPosition(position);
    };

    if (isLoading) return (
        <div className="w-full h-full flex items-center justify-center bg-gray-50/50 backdrop-blur-sm">
            <div className="flex flex-col items-center gap-4">
                <div className="w-12 h-12 border-4 border-blue-500 border-t-transparent rounded-full animate-spin" />
                <span className="text-sm font-medium text-gray-500">Aligning Molecular Structures...</span>
            </div>
        </div>
    );

    if (!dashboard) return null;

    const propertyMetrics = [
        { key: 'molecular_weight' as const, id: 'MW', label: 'Mol Weight', unit: 'g' },
        { key: 'logp' as const, id: 'logp', label: 'LogP', unit: '' },
        { key: 'tpsa' as const, id: 'tpsa', label: 'TPSA', unit: 'A2' },
        { key: 'qed' as const, id: 'qed', label: 'QED', unit: '' }
    ];

    return (
        <div className="flex flex-col lg:flex-row w-full h-full gap-4 p-4 animate-in fade-in duration-700">
            <div className="flex-[2] flex flex-col gap-4">
                <div className="grid grid-cols-2 gap-4 flex-1">
                    <div className="relative rounded-[2rem] overflow-hidden border border-gray-100 bg-white shadow-sm group">
                        <div className="absolute top-6 left-6 z-10 flex items-center gap-2 px-3 py-1.5 bg-white/80 backdrop-blur-md rounded-full border border-gray-100 shadow-sm">
                            <div className="w-2 h-2 rounded-full bg-gray-400" />
                            <span className="text-[10px] font-black uppercase tracking-widest text-gray-500">Baseline</span>
                        </div>
                        <Studio3DScene
                            mode={mode}
                            molecule={dashboard.base_version.json_graph}
                            diffData={dashboard.diff.baseline}
                            onCameraChange={handleCameraChange}
                            cameraTarget={cameraTarget}
                            cameraPosition={cameraPosition}
                        />
                    </div>

                    <div className="relative rounded-[2rem] overflow-hidden border-2 border-blue-100 bg-white shadow-xl shadow-blue-500/5 group">
                        <div className="absolute top-6 left-6 z-10 flex items-center gap-2 px-3 py-1.5 bg-blue-600 rounded-full shadow-lg shadow-blue-500/20">
                            <div className="w-2 h-2 rounded-full bg-white animate-pulse" />
                            <span className="text-[10px] font-black uppercase tracking-widest text-white">Optimized Proposal</span>
                        </div>
                        <Studio3DScene
                            mode={mode}
                            molecule={dashboard.opt_version.json_graph}
                            diffData={dashboard.diff.proposal}
                            onCameraChange={handleCameraChange}
                            cameraTarget={cameraTarget}
                            cameraPosition={cameraPosition}
                        />
                    </div>
                </div>

                <div className="h-48 bg-white/50 backdrop-blur-md rounded-[2rem] border border-gray-100 p-6 overflow-hidden flex gap-8">
                    <div className="flex-1">
                        <h4 className="text-[10px] font-black uppercase tracking-widest text-gray-400 mb-4 flex items-center gap-2">
                            <ArrowRight className="w-3 h-3" /> Property Deltas
                        </h4>
                        <table className="w-full text-left">
                            <thead>
                                <tr className="text-[10px] font-bold text-gray-400 uppercase tracking-tighter">
                                    <th className="pb-2">Metric</th>
                                    <th className="pb-2">Base</th>
                                    <th className="pb-2">Prop</th>
                                    <th className="pb-2">Delta</th>
                                </tr>
                            </thead>
                            <tbody className="text-xs font-medium">
                                {propertyMetrics.map((metric) => {
                                    const baseValue = getNumericProperty(dashboard.base_version.properties, metric.key);
                                    const proposalValue = getNumericProperty(dashboard.opt_version.properties, metric.key);
                                    const delta = proposalValue - baseValue;

                                    return (
                                        <tr key={metric.id} className="border-t border-gray-50">
                                            <td className="py-2 text-gray-500">{metric.label}</td>
                                            <td className="py-2">{baseValue.toFixed(2)}</td>
                                            <td className="py-2">{proposalValue.toFixed(2)}</td>
                                            <td className={`py-2 ${delta > 0 ? 'text-blue-600' : 'text-emerald-600'}`}>
                                                {delta > 0 ? '+' : ''}{delta.toFixed(2)}{metric.unit}
                                            </td>
                                        </tr>
                                    );
                                })}
                            </tbody>
                        </table>
                    </div>

                    <div className="w-px bg-gray-100" />

                    <div className="flex-1">
                        <h4 className="text-[10px] font-black uppercase tracking-widest text-gray-400 mb-2 flex items-center gap-2">
                            <Zap className="w-3 h-3 text-amber-500" /> Profile Shift
                        </h4>
                        <div className="h-32">
                            <RadarChart baseline={dashboard.radar} />
                        </div>
                    </div>
                </div>
            </div>

            <div className="flex-1 flex flex-col gap-4">
                <div className="bg-white rounded-[2rem] border border-gray-100 p-6 shadow-sm flex flex-col gap-4 flex-1">
                    <div>
                        <h3 className="text-sm font-black uppercase tracking-widest text-gray-800 mb-4">Quality Control</h3>
                        <div className="space-y-3">
                            {dashboard.optimization_issues.map((issue: AnalysisIssue) => (
                                <div key={`${issue.title}-${issue.description}`} className="p-3 bg-red-50 rounded-xl border border-red-100 flex gap-3">
                                    <AlertTriangle className="w-4 h-4 text-red-500 shrink-0" />
                                    <div>
                                        <div className="text-[10px] font-bold text-red-800 uppercase">{issue.title}</div>
                                        <div className="text-[10px] text-red-600/80 leading-tight">{issue.description}</div>
                                    </div>
                                </div>
                            ))}
                            {dashboard.optimization_issues.length === 0 && (
                                <div className="p-3 bg-emerald-50 rounded-xl border border-emerald-100 flex gap-3">
                                    <CheckCircle className="w-4 h-4 text-emerald-500 shrink-0" />
                                    <div className="text-[10px] font-bold text-emerald-800 uppercase">Structural Integrity Verified</div>
                                </div>
                            )}
                        </div>
                    </div>

                    <div className="mt-auto pt-6 border-t border-gray-100">
                        <button className="w-full py-4 bg-blue-600 text-white rounded-2xl font-bold text-xs uppercase tracking-widest shadow-lg shadow-blue-500/30 hover:bg-blue-700 transition-all active:scale-[0.98]">
                            Accept Optimization
                        </button>
                        <button className="w-full py-4 bg-transparent text-gray-400 rounded-2xl font-bold text-xs uppercase tracking-widest mt-2 hover:bg-gray-50 transition-all">
                            Reject Proposal
                        </button>
                    </div>
                </div>
            </div>
        </div>
    );
}
