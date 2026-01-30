import { useState } from 'react';
import { useStudioStore } from '../../store/studioStore';
import Studio3DScene from './Studio3DScene';
import * as THREE from 'three';
import { BoxSelect, Maximize2, Layers } from 'lucide-react';

export default function StudioMainCanvas() {
    const { status, dashboard } = useStudioStore();
    const [cameraState, setCameraState] = useState({
        target: new THREE.Vector3(0, 0, 0),
        position: new THREE.Vector3(0, 0, 8)
    });

    const handleCameraChange = (target: THREE.Vector3, position: THREE.Vector3) => {
        setCameraState({ target, position });
    };

    if (!dashboard) return (
        <div className="flex-1 flex items-center justify-center bg-[#F9FAFB]">
            <div className="animate-pulse flex flex-col items-center gap-3">
                <div className="w-12 h-12 rounded-full border-4 border-blue-100 border-t-blue-500" />
                <span className="text-[10px] font-black text-gray-300 uppercase tracking-widest">Waking Kernel...</span>
            </div>
        </div>
    );

    const isInDiffMode = status === 'PROPOSED' || status === 'COMMITTED';

    return (
        <div className="flex-1 flex divide-x divide-[#E5E7EB] bg-white relative">
            {/* Viewport 1: Baseline */}
            <div className="flex-1 relative group overflow-hidden border-t border-[#E5E7EB]">
                {/* Truth Pane Header */}
                <div className="absolute top-0 left-0 right-0 h-8 bg-slate-50/80 backdrop-blur-sm border-b border-[#E5E7EB] z-10 flex items-center justify-between px-3">
                    <div className="flex items-center gap-2">
                        <div className="w-1.5 h-1.5 rounded-full bg-blue-500 shadow-[0_0_8px_rgba(59,130,246,0.5)]" />
                        <span className="text-[9px] font-black text-gray-900 uppercase tracking-widest">Structural Truth</span>
                        <span className="text-[8px] font-medium text-gray-400 font-mono">v.BASELINE</span>
                    </div>
                    <div className="flex gap-2">
                        <button className="text-gray-300 hover:text-gray-500 transition-colors"><BoxSelect size={12} /></button>
                        <button className="text-gray-300 hover:text-gray-500 transition-colors"><Layers size={12} /></button>
                    </div>
                </div>

                <div className="w-full h-full pt-8">
                    <Studio3DScene
                        mode="optimize"
                        molecule={dashboard.baseline.graph}
                        cameraTarget={cameraState.target}
                        cameraPosition={cameraState.position}
                        onCameraChange={handleCameraChange}
                        diffData={dashboard.diff as any}
                    />
                </div>
            </div>

            {/* Viewport 2: Proposal (Conditional) */}
            {isInDiffMode && (
                <div className="flex-1 relative group bg-blue-50/10 border-t border-[#E5E7EB] animate-in slide-in-from-right duration-500">
                    {/* Potential Pane Header */}
                    <div className="absolute top-0 left-0 right-0 h-8 bg-blue-600 z-10 flex items-center justify-between px-3 shadow-md">
                        <div className="flex items-center gap-2">
                            <div className="w-1.5 h-1.5 rounded-full bg-white shadow-[0_0_8px_rgba(255,255,255,0.8)] animate-pulse" />
                            <span className="text-[9px] font-black text-white uppercase tracking-widest leading-none">AI Optimized Proposal</span>
                            <span className="text-[8px] font-medium text-blue-200 font-mono">v.EPHEMERAL</span>
                        </div>
                        <button className="text-blue-200 hover:text-white transition-colors"><Maximize2 size={12} /></button>
                    </div>

                    <div className="w-full h-full pt-8">
                        <Studio3DScene
                            mode="optimize"
                            molecule={dashboard.proposal?.graph || null}
                            cameraTarget={cameraState.target}
                            cameraPosition={cameraState.position}
                            onCameraChange={handleCameraChange}
                            diffData={dashboard.diff as any}
                        />
                    </div>
                </div>
            )}
        </div>
    );
}
