
import { useStudioStore } from '../../store/studioStore';
import { Activity, ShieldCheck } from 'lucide-react';

export default function StudioTopBar() {
    const { status, baselineVersionId, proposalVersionId } = useStudioStore();

    return (
        <header className="h-14 bg-white border-b border-[#E5E7EB] flex items-center justify-between px-6 shrink-0">
            <div className="flex items-center gap-6">
                <div className="flex items-center gap-2">
                    <Activity size={18} className="text-blue-600" />
                    <h1 className="text-sm font-black uppercase tracking-widest text-black">Studio</h1>
                </div>

                <div className="h-4 w-[1px] bg-[#E5E7EB]" />

                <div className="flex items-center gap-4">
                    <div className="flex flex-col">
                        <span className="text-[8px] font-black text-gray-400 uppercase tracking-tighter">Baseline</span>
                        <span className="text-[10px] font-mono text-gray-600">{baselineVersionId?.slice(0, 8) || 'N/A'}</span>
                    </div>
                    {proposalVersionId && (
                        <>
                            <div className="h-2 w-2 rounded-full bg-gray-300" />
                            <div className="flex flex-col">
                                <span className="text-[8px] font-black text-gray-400 uppercase tracking-tighter">Proposal</span>
                                <span className="text-[10px] font-mono text-blue-600 font-bold">{proposalVersionId.slice(0, 8)}</span>
                            </div>
                        </>
                    )}
                </div>
            </div>

            <div className="flex items-center gap-3">
                <div className={`flex items-center gap-2 px-3 py-1 rounded-full border ${status === 'READY' ? 'bg-green-50 border-green-200 text-green-700' :
                    status === 'PROPOSED' ? 'bg-blue-50 border-blue-200 text-blue-700' :
                        status === 'OPTIMIZING' ? 'bg-amber-50 border-amber-200 text-amber-700' :
                            'bg-gray-50 border-gray-200 text-gray-600'
                    }`}>
                    <div className={`w-1.5 h-1.5 rounded-full ${status === 'READY' ? 'bg-green-500' :
                        status === 'PROPOSED' ? 'bg-blue-500' :
                            status === 'OPTIMIZING' ? 'bg-amber-500 animate-pulse' :
                                'bg-gray-400'
                        }`} />
                    <span className="text-[9px] font-black uppercase tracking-widest">{status}</span>
                </div>

                <div className="flex items-center gap-1.5 px-3 py-1 rounded-md bg-[#F9FAFB] border border-[#E5E7EB]">
                    <ShieldCheck size={12} className="text-gray-400" />
                    <span className="text-[9px] font-bold text-gray-500">DETERMINISTIC KERNEL</span>
                </div>
            </div>
        </header>
    );
}
