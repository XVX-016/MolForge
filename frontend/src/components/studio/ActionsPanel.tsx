
import { useState, useRef, useEffect, useMemo } from 'react';
import { useStudioStore } from '../../store/studioStore';
import { useHistoryStore } from '../../store/historyStore';
import { useStudioMode } from '../../lib/studio/hooks';
import { Send, Clock, Activity, FileJson, FileType, RotateCcw, Sparkles, Plus } from 'lucide-react';
import { useStudioGateway } from '../../lib/ai/studioGateway';
import Card from '../ui/Card';
import { exportToJSON, exportToSMILES, downloadFile } from '../../lib/io';
import { calculateMolecularWeight, estimateLogP } from '../../lib/chemistry';

export default function ActionsPanel() {
    const [inputValue, setInputValue] = useState('');
    const [activeTab, setActiveTab] = useState<'actions' | 'logs' | 'compare'>('actions');

    const bottomRef = useRef<HTMLDivElement>(null);
    const { messages, addMessage, mode } = useStudioStore();
    const { present: molecule, past, logs } = useHistoryStore();
    const { sendCommand, isProcessing } = useStudioGateway();
    const { modeColor, canSimulate } = useStudioMode();

    useEffect(() => {
        bottomRef.current?.scrollIntoView({ behavior: 'smooth' });
    }, [messages, logs]);

    const handleActionClick = async (cmd: string) => {
        if (!molecule) return;
        addMessage({
            id: Date.now().toString(),
            role: 'user',
            content: cmd,
            timestamp: new Date()
        });
        await sendCommand(cmd, mode, molecule);
    };

    const handleChatSubmit = async (e: React.FormEvent) => {
        e.preventDefault();
        if (!inputValue.trim() || !molecule) return;
        handleActionClick(inputValue);
        setInputValue('');
    };

    // --- Comparison Logic ---
    const comparisonData = useMemo(() => {
        if (!molecule) return null;
        const initial = past.length > 0 ? past[0] : molecule;
        const distinct = past.length > 0;
        return {
            distinct,
            currentStats: { weight: calculateMolecularWeight(molecule), logP: estimateLogP(molecule) },
            initialStats: { weight: calculateMolecularWeight(initial), logP: estimateLogP(initial) }
        };
    }, [molecule, past]);

    return (
        <Card className="flex flex-col h-full bg-white border-lightGrey shadow-sm overflow-hidden">
            {/* Header Tabs */}
            <div className="flex items-center gap-1 p-2 border-b border-lightGrey bg-gray-50/50">
                <button onClick={() => setActiveTab('actions')} className={`px-3 py-1.5 rounded-lg text-xs font-bold transition-all flex items-center gap-2 ${activeTab === 'actions' ? 'bg-white text-black shadow-sm border border-gray-200' : 'text-midGrey hover:bg-gray-100'}`}>
                    <Activity size={14} /> Actions
                </button>
                <button onClick={() => setActiveTab('logs')} className={`px-3 py-1.5 rounded-lg text-xs font-bold transition-all flex items-center gap-2 ${activeTab === 'logs' ? 'bg-white text-black shadow-sm border border-gray-200' : 'text-midGrey hover:bg-gray-100'}`}>
                    <Clock size={14} /> Log
                </button>
                <button onClick={() => setActiveTab('compare')} className={`px-3 py-1.5 rounded-lg text-xs font-bold transition-all flex items-center gap-2 ${activeTab === 'compare' ? 'bg-white text-black shadow-sm border border-gray-200' : 'text-midGrey hover:bg-gray-100'}`}>
                    <RotateCcw size={14} /> Compare
                </button>
                <div className="flex-1" />
                <div className="flex items-center gap-1">
                    <button onClick={() => downloadFile(`mol_${Date.now()}.smi`, exportToSMILES(molecule!), 'text/plain')} className="p-1.5 rounded-md hover:bg-gray-200 text-midGrey hover:text-black transition-colors" title="Export SMILES"><FileType size={14} /></button>
                    <button onClick={() => downloadFile(`mol_${Date.now()}.json`, exportToJSON(molecule!), 'application/json')} className="p-1.5 rounded-md hover:bg-gray-200 text-midGrey hover:text-black transition-colors" title="Export JSON"><FileJson size={14} /></button>
                </div>
            </div>

            {/* CONTENT AREA */}
            <div className="flex-1 overflow-y-auto p-5 space-y-6 bg-white relative">

                {activeTab === 'actions' && (
                    <div className="animate-in fade-in slide-in-from-right-2 duration-300 space-y-6">
                        {/* MODE HEADER */}
                        <div>
                            <div className="flex items-center gap-2 mb-1">
                                <h2 className="text-sm font-black uppercase tracking-[0.15em]" style={{ color: modeColor }}>
                                    {mode} Mode
                                </h2>
                            </div>
                            <p className="text-[11px] text-midGrey leading-relaxed">
                                {mode === 'design' ? 'Structural Authoring Interface. Direct manipulation of atoms and bonds is active.' :
                                    mode === 'optimize' ? 'Optimization Decision Hub. AI-driven suggestions for property improvement.' :
                                        'Simulation Playback Environment. Dynamics analysis and energy monitoring.'}
                            </p>
                        </div>

                        {/* PRIMARY ACTIONS */}
                        <div className="space-y-3">
                            <h3 className="text-[10px] font-bold text-gray-400 uppercase tracking-widest">Primary Actions</h3>
                            <div className="grid grid-cols-1 gap-2">
                                {mode === 'design' && (
                                    <>
                                        <button onClick={() => handleActionClick("Add benzene ring")} className="flex items-center justify-between px-4 py-3 bg-gray-50 hover:bg-black hover:text-white border border-lightGrey rounded-xl transition-all group text-xs font-bold">
                                            <span>Add Benzene Ring</span>
                                            <Plus size={14} className="text-gray-300 group-hover:text-white" />
                                        </button>
                                        <button onClick={() => handleActionClick("Validate structure")} className="flex items-center justify-between px-4 py-3 bg-gray-50 hover:bg-black hover:text-white border border-lightGrey rounded-xl transition-all group text-xs font-bold">
                                            <span>Check Valence Health</span>
                                            <Activity size={14} className="text-gray-300 group-hover:text-white" />
                                        </button>
                                    </>
                                )}
                                {mode === 'optimize' && (
                                    <button className="flex items-center justify-between px-4 py-3 bg-purple-50 text-purple-700 hover:bg-purple-600 hover:text-white border border-purple-100 rounded-xl transition-all group text-xs font-bold">
                                        <span>Scan for Improvements</span>
                                        <Sparkles size={14} className="text-purple-300 group-hover:text-white" />
                                    </button>
                                )}
                                {mode === 'simulate' && (
                                    <div className="p-4 bg-orange-50 rounded-xl border border-orange-100 text-[11px] text-orange-800 font-medium leading-relaxed italic">
                                        "Simulation mode is currently set to observational playback. Interaction is restricted to the timeline."
                                    </div>
                                )}
                            </div>
                        </div>

                        {/* SECONDARY/AI ACTIONS */}
                        {!canSimulate && (
                            <div className="space-y-3 pt-2">
                                <h3 className="text-[10px] font-bold text-gray-400 uppercase tracking-widest">AI Assist (Optional)</h3>
                                <form onSubmit={handleChatSubmit} className="relative group">
                                    <input
                                        type="text"
                                        value={inputValue}
                                        onChange={(e) => setInputValue(e.target.value)}
                                        placeholder="Describe a targeted change..."
                                        className="w-full bg-offwhite border border-lightGrey rounded-xl pl-4 pr-12 py-3 text-[11px] focus:outline-none focus:ring-2 focus:ring-black/5 focus:border-black/20 transition-all font-medium"
                                        disabled={isProcessing}
                                    />
                                    <button type="submit" disabled={!inputValue.trim() || isProcessing} className="absolute right-2 top-1/2 -translate-y-1/2 p-1.5 bg-black text-white rounded-lg hover:bg-gray-800 disabled:opacity-20 transition-all">
                                        {isProcessing ? <div className="w-4 h-4 border-2 border-white/30 border-t-white rounded-full animate-spin" /> : <Send size={14} />}
                                    </button>
                                </form>
                            </div>
                        )}
                    </div>
                )}

                {/* TAB: LOGS */}
                {activeTab === 'logs' && (
                    <div className="space-y-4 animate-in fade-in duration-300">
                        <div className="flex items-center justify-between">
                            <h3 className="text-[10px] font-bold text-gray-400 uppercase tracking-widest">Chronological Action Log</h3>
                            <button className="text-[10px] text-blue-600 font-bold hover:underline">Clear History</button>
                        </div>
                        <div className="space-y-3">
                            {logs.length === 0 && <p className="text-xs text-center text-midGrey italic mt-10">No actions recorded yet.</p>}
                            {logs.map((log) => (
                                <div key={log.id} className="group relative pl-4 border-l-2 border-lightGrey hover:border-black transition-colors py-1">
                                    <div className="flex items-center justify-between mb-0.5">
                                        <span className="text-[11px] font-bold text-black">{log.description}</span>
                                        <span className="text-[9px] font-mono text-midGrey">{new Date(log.timestamp).toLocaleTimeString([], { hour: '2-digit', minute: '2-digit', second: '2-digit' })}</span>
                                    </div>
                                    <div className="text-[9px] text-midGrey uppercase tracking-tighter flex items-center gap-2">
                                        <span className="bg-gray-100 px-1 rounded">{log.source}</span>
                                        <span>/</span>
                                        <span style={{ color: log.mode === 'design' ? '#22c55e' : log.mode === 'optimize' ? '#a855f7' : '#f97316' }}>{log.mode}</span>
                                    </div>
                                </div>
                            ))}
                        </div>
                        <div ref={bottomRef} />
                    </div>
                )}

                {/* TAB: COMPARE */}
                {activeTab === 'compare' && comparisonData && (
                    <div className="animate-in fade-in duration-300 space-y-6">
                        <div className="bg-gray-50 p-5 rounded-2xl border border-lightGrey/50 text-center shadow-sm">
                            <p className="text-[10px] text-midGrey uppercase tracking-[0.2em] font-black mb-1">State Confidence</p>
                            {comparisonData.distinct ? (
                                <span className="text-sm font-bold text-blue-600">MODIFIED (± {past.length} DELTAS)</span>
                            ) : (
                                <span className="text-sm font-bold text-gray-400">INITIAL STATE</span>
                            )}
                        </div>

                        <div className="grid grid-cols-3 gap-y-6 items-center">
                            <div className="text-[9px] font-black text-midGrey uppercase tracking-widest">Metric</div>
                            <div className="text-[9px] font-black text-midGrey uppercase tracking-widest text-center">Initial</div>
                            <div className="text-[9px] font-black text-midGrey uppercase tracking-widest text-center">Current</div>

                            <div className="text-xs font-bold text-black">Molecular Weight</div>
                            <div className="text-xs font-mono text-midGrey text-center">{comparisonData.initialStats.weight}</div>
                            <div className={`text-xs font-mono font-bold text-center ${comparisonData.currentStats.weight !== comparisonData.initialStats.weight ? 'text-blue-600' : 'text-black'}`}>{comparisonData.currentStats.weight}</div>

                            <div className="text-xs font-bold text-black">Estimated LogP</div>
                            <div className="text-xs font-mono text-midGrey text-center">{comparisonData.initialStats.logP}</div>
                            <div className={`text-xs font-mono font-bold text-center ${comparisonData.currentStats.logP !== comparisonData.initialStats.logP ? 'text-blue-600' : 'text-black'}`}>{comparisonData.currentStats.logP}</div>
                        </div>

                        {comparisonData.distinct && (
                            <div className="p-4 bg-blue-50/50 border border-blue-100 rounded-xl text-[11px] text-blue-800 leading-relaxed shadow-sm">
                                <strong className="block mb-1 uppercase tracking-wider text-[9px]">Divergence Summary</strong>
                                Structure has drifted by <b>{(comparisonData.currentStats.weight - comparisonData.initialStats.weight).toFixed(2)} Da</b> and <b>{(comparisonData.currentStats.logP - comparisonData.initialStats.logP).toFixed(2)} LogP</b> from root.
                            </div>
                        )}
                    </div>
                )}
            </div>
        </Card>
    );
}
