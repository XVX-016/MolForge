
import { useState, useRef, useEffect, useMemo } from 'react';
import { useStudioStore } from '../../store/studioStore';
import { useHistoryStore } from '../../store/historyStore';
import { Send, Terminal, Clock, Activity, FileJson, FileType } from 'lucide-react';
import { useStudioGateway } from '../../lib/ai/studioGateway';
import Card from '../ui/Card';
import { exportToJSON, exportToSMILES, downloadFile } from '../../lib/io';
import { calculateMolecularWeight, estimateLogP } from '../../lib/chemistry';

export default function CommandPanel() {
    const [inputValue, setInputValue] = useState('');
    const [activeTab, setActiveTab] = useState<'command' | 'logs' | 'compare'>('command');

    const bottomRef = useRef<HTMLDivElement>(null);
    const { messages, addMessage, mode } = useStudioStore();
    const { present: molecule, past, logs } = useHistoryStore();
    const { sendCommand, isProcessing } = useStudioGateway();

    // --- Auto-scroll ---
    useEffect(() => {
        bottomRef.current?.scrollIntoView({ behavior: 'smooth' });
    }, [messages, logs]);

    const handleSubmit = async (e: React.FormEvent) => {
        e.preventDefault();
        if (!inputValue.trim() || !molecule) return;

        // Optimistic UI update
        addMessage({
            id: Date.now().toString(),
            role: 'user',
            content: inputValue,
            timestamp: Date.now()
        });

        const command = inputValue;
        setInputValue(''); // Clear input immediately

        await sendCommand(command, mode, molecule);
    };

    // --- Export Actions ---
    const handleExportJSON = () => {
        if (!molecule) return;
        const json = exportToJSON(molecule);
        downloadFile(`molecule_${Date.now()}.json`, json, 'application/json');
    };

    const handleExportSMILES = () => {
        if (!molecule) return;
        const smiles = exportToSMILES(molecule);
        downloadFile(`molecule_${Date.now()}.smi`, smiles, 'text/plain');
    };

    // --- Comparison Logic ---
    const comparisonData = useMemo(() => {
        if (!molecule) return null;
        const initial = past.length > 0 ? past[0] : molecule;

        const distinct = past.length > 0; // Truly different history exists

        const currentStats = {
            weight: calculateMolecularWeight(molecule),
            logP: estimateLogP(molecule)
        };

        const initialStats = {
            weight: calculateMolecularWeight(initial),
            logP: estimateLogP(initial)
        };

        return { distinct, currentStats, initialStats };
    }, [molecule, past]);

    return (
        <Card className="flex flex-col h-full bg-white border-lightGrey shadow-sm overflow-hidden">

            {/* Header Tabs */}
            <div className="flex items-center gap-1 p-2 border-b border-lightGrey bg-gray-50/50">
                <button
                    onClick={() => setActiveTab('command')}
                    className={`px-3 py-1.5 rounded-lg text-xs font-bold transition-all flex items-center gap-2 ${activeTab === 'command'
                            ? 'bg-white text-black shadow-sm border border-gray-200'
                            : 'text-midGrey hover:bg-gray-100'
                        }`}
                >
                    <Terminal size={14} /> Command
                </button>
                <button
                    onClick={() => setActiveTab('logs')}
                    className={`px-3 py-1.5 rounded-lg text-xs font-bold transition-all flex items-center gap-2 ${activeTab === 'logs'
                            ? 'bg-white text-black shadow-sm border border-gray-200'
                            : 'text-midGrey hover:bg-gray-100'
                        }`}
                >
                    <Clock size={14} /> Log
                </button>
                <button
                    onClick={() => setActiveTab('compare')}
                    className={`px-3 py-1.5 rounded-lg text-xs font-bold transition-all flex items-center gap-2 ${activeTab === 'compare'
                            ? 'bg-white text-black shadow-sm border border-gray-200'
                            : 'text-midGrey hover:bg-gray-100'
                        }`}
                >
                    <Activity size={14} /> Compare
                </button>

                <div className="flex-1" />

                {/* Mini Export Dropdown/Buttons */}
                <div className="flex items-center gap-1">
                    <button
                        onClick={handleExportSMILES}
                        className="p-1.5 rounded-md hover:bg-gray-200 text-midGrey hover:text-black transition-colors"
                        title="Export SMILES"
                    >
                        <FileType size={14} />
                    </button>
                    <button
                        onClick={handleExportJSON}
                        className="p-1.5 rounded-md hover:bg-gray-200 text-midGrey hover:text-black transition-colors"
                        title="Export JSON"
                    >
                        <FileJson size={14} />
                    </button>
                </div>
            </div>

            {/* CONTENT AREA */}
            <div className="flex-1 overflow-y-auto p-4 space-y-4 bg-white relative">

                {/* TAB: COMMAND */}
                {activeTab === 'command' && (
                    <>
                        {messages.length === 0 ? (
                            <div className="h-full flex flex-col items-center justify-center text-center opacity-40">
                                <Terminal size={32} className="mb-2" />
                                <p className="text-sm font-medium">Ready for commands</p>
                                <p className="text-xs">Try "Optimize for solubility"</p>
                            </div>
                        ) : (
                            messages.map((msg) => (
                                <div key={msg.id} className={`flex ${msg.role === 'user' ? 'justify-end' : 'justify-start'}`}>
                                    <div className={`
                                max-w-[85%] rounded-2xl px-4 py-2.5 text-sm leading-relaxed shadow-sm
                                ${msg.role === 'user'
                                            ? 'bg-black text-white rounded-br-none'
                                            : 'bg-offwhite border border-lightGrey text-darkGrey rounded-bl-none'}
                            `}>
                                        {msg.content}
                                    </div>
                                </div>
                            ))
                        )}
                        <div ref={bottomRef} />
                    </>
                )}

                {/* TAB: LOGS */}
                {activeTab === 'logs' && (
                    <div className="space-y-3">
                        {logs.length === 0 && (
                            <p className="text-xs text-center text-midGrey italic mt-10">No actions recorded yet.</p>
                        )}
                        {logs.map((log) => (
                            <div key={log.id} className="flex gap-3 text-xs border-b border-lightGrey/50 pb-2 last:border-0">
                                <div className="font-mono text-midGrey min-w-[50px]">
                                    {new Date(log.timestamp).toLocaleTimeString([], { hour: '2-digit', minute: '2-digit', second: '2-digit' })}
                                </div>
                                <div className="flex-1">
                                    <div className="font-medium text-black">{log.description}</div>
                                    <div className="text-[10px] text-gray-400 capitalize mt-0.5 flex items-center gap-1">
                                        {log.source} • <span className="uppercase tracking-wider">{log.mode}</span>
                                    </div>
                                </div>
                            </div>
                        ))}
                        <div ref={bottomRef} />
                    </div>
                )}

                {/* TAB: COMPARE */}
                {activeTab === 'compare' && comparisonData && (
                    <div className="flex flex-col gap-6">
                        <div className="bg-gray-50 p-4 rounded-xl border border-lightGrey text-center">
                            <p className="text-xs text-midGrey uppercase tracking-wider font-bold mb-1">Status</p>
                            {comparisonData.distinct ? (
                                <span className="text-sm font-bold text-blue-600">Modified ({past.length} Steps)</span>
                            ) : (
                                <span className="text-sm font-bold text-gray-400">Unchanged</span>
                            )}
                        </div>

                        {/* Diff Table */}
                        <div className="grid grid-cols-3 gap-y-4 items-center text-center">
                            {/* Headings */}
                            <div className="text-[10px] font-bold text-midGrey uppercase">Metric</div>
                            <div className="text-[10px] font-bold text-midGrey uppercase">Initial</div>
                            <div className="text-[10px] font-bold text-midGrey uppercase">Current</div>

                            {/* Weight */}
                            <div className="text-xs font-bold text-black text-left pl-2">Weight</div>
                            <div className="text-xs font-mono text-midGrey">{comparisonData.initialStats.weight}</div>
                            <div className={`text-xs font-mono font-bold ${comparisonData.currentStats.weight !== comparisonData.initialStats.weight ? 'text-blue-600' : 'text-black'}`}>
                                {comparisonData.currentStats.weight}
                            </div>

                            {/* LogP */}
                            <div className="text-xs font-bold text-black text-left pl-2">Est. LogP</div>
                            <div className="text-xs font-mono text-midGrey">{comparisonData.initialStats.logP}</div>
                            <div className={`text-xs font-mono font-bold ${comparisonData.currentStats.logP !== comparisonData.initialStats.logP ? 'text-blue-600' : 'text-black'}`}>
                                {comparisonData.currentStats.logP}
                            </div>
                        </div>

                        {comparisonData.distinct && (
                            <div className="mt-4 p-3 bg-blue-50 border border-blue-100 rounded-lg text-xs text-blue-800 leading-relaxed">
                                <strong className="block mb-1">Change Summary</strong>
                                Weight shifted by <b>{(comparisonData.currentStats.weight - comparisonData.initialStats.weight).toFixed(2)}</b>.<br />
                                LogP shifted by <b>{(comparisonData.currentStats.logP - comparisonData.initialStats.logP).toFixed(2)}</b>.
                            </div>
                        )}
                    </div>
                )}

            </div>

            {/* FOOTER: INPUT (Only for Command Tab) */}
            {activeTab === 'command' && (
                <div className="p-3 bg-white border-t border-lightGrey shrink-0 z-20">
                    <form onSubmit={handleSubmit} className="relative group">
                        <input
                            type="text"
                            value={inputValue}
                            onChange={(e) => setInputValue(e.target.value)}
                            placeholder={mode === 'design' ? "Describe a change..." : "Enter command..."}
                            className="w-full bg-offwhite border border-lightGrey rounded-xl pl-4 pr-12 py-3 text-sm focus:outline-none focus:ring-2 focus:ring-black/5 focus:border-black/20 transition-all placeholder:text-midGrey"
                            disabled={isProcessing}
                        />
                        <button
                            type="submit"
                            disabled={!inputValue.trim() || isProcessing}
                            className="absolute right-2 top-1/2 -translate-y-1/2 p-1.5 bg-black text-white rounded-lg hover:bg-gray-800 disabled:opacity-50 disabled:hover:bg-black transition-all"
                        >
                            {isProcessing ? (
                                <div className="w-4 h-4 border-2 border-white/30 border-t-white rounded-full animate-spin" />
                            ) : (
                                <Send size={16} />
                            )}
                        </button>
                    </form>
                </div>
            )}
        </Card>
    );
}
