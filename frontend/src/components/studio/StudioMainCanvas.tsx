
import { useState } from 'react';
import { useStudioStore } from '../../store/studioStore';
import Studio3DScene from './Studio3DScene';
import * as THREE from 'three';

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
        <div className="flex-1 flex divide-x divide-[#E5E7EB] bg-white">
            {/* Viewport 1: Baseline */}
            <div className="flex-1 relative group">
                <div className="absolute top-4 left-4 z-10">
                    <span className="px-2 py-1 bg-white/80 backdrop-blur-sm border border-[#E5E7EB] rounded text-[8px] font-black text-gray-400 uppercase tracking-widest">
                        Baseline
                    </span>
                </div>
                <Studio3DScene
                    mode="optimize"
                    molecule={dashboard.baseline.graph}
                    cameraTarget={cameraState.target}
                    cameraPosition={cameraState.position}
                    onCameraChange={handleCameraChange}
                    diffData={dashboard.diff as any}
                />
            </div>

            {/* Viewport 2: Proposal (Conditional) */}
            {isInDiffMode && (
                <div className="flex-1 relative group bg-blue-50/10">
                    <div className="absolute top-4 left-4 z-10">
                        <span className="px-2 py-1 bg-blue-600 text-white rounded text-[8px] font-black uppercase tracking-widest shadow-lg shadow-blue-200">
                            Proposal
                        </span>
                    </div>
                    <Studio3DScene
                        mode="optimize"
                        molecule={dashboard.proposal?.graph || null}
                        cameraTarget={cameraState.target}
                        cameraPosition={cameraState.position}
                        onCameraChange={handleCameraChange}
                        diffData={dashboard.diff as any}
                    />
                </div>
            )}
        </div>
    );
}
