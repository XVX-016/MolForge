import { useState, useRef, useEffect } from 'react';
import { MessageSquare, Send, Sparkles, Loader2 } from 'lucide-react';
import { useStudioStore } from '../../store/studioStore';
import { useHistoryStore } from '../../store/historyStore';
import { motion, AnimatePresence } from 'framer-motion';

type PanelMode = 'CHAT' | 'LOG';

export default function AIPanel() {
    const [panelMode, setPanelMode] = useState<PanelMode>('CHAT');
    const [command, setCommand] = useState('');
    const { isCommandRunning, runCommand } = useStudioStore();
    const { logs } = useHistoryStore();
    const scrollRef = useRef<HTMLDivElement>(null);

    // Auto-scroll logs/messages
    useEffect(() => {
        if (scrollRef.current) {
            scrollRef.current.scrollTop = scrollRef.current.scrollHeight;
        }
    }, [logs, panelMode]);

    const handleSubmit = async (e: React.FormEvent) => {
        e.preventDefault();
        if (!command.trim() || isCommandRunning) return;

        const currentCommand = command;
        setCommand('');
        await runCommand(currentCommand);
    };

    const suggestions = [
        { label: 'Create Benzene', icon: '⬡' },
        { label: 'Add Nitro groups', icon: 'NO₂' },
        { label: 'Optimize Geometry', icon: '✨' },
    ];

    return (
        <div className="flex flex-col h-full bg-white">
            {/* Panel Header */}
            <div className="flex items-center justify-between p-4 border-b border-[#E5E7EB] bg-[#F9FAFB]">
                <div className="flex items-center gap-2">
                    <Sparkles size={14} className="text-blue-600" />
                    <span className="text-[10px] font-black uppercase tracking-[0.2em] text-[#6B7280]">
                        Architect <span className="text-blue-600/80">v1.5</span>
                    </span>
                </div>

                <div className="flex bg-[#F3F4F6] rounded-md p-0.5 border border-[#E5E7EB]">
                    <button
                        onClick={() => setPanelMode('CHAT')}
                        className={`px-3 py-1 rounded-sm text-[9px] font-bold transition-all ${panelMode === 'CHAT' ? 'bg-white text-black shadow-sm' : 'text-gray-500 hover:text-gray-700'}`}
                    >
                        CHAT
                    </button>
                    <button
                        onClick={() => setPanelMode('LOG')}
                        className={`px-3 py-1 rounded-sm text-[9px] font-bold transition-all ${panelMode === 'LOG' ? 'bg-white text-black shadow-sm' : 'text-gray-500 hover:text-gray-700'}`}
                    >
                        LOGS
                    </button>
                </div>
            </div>

            {/* Panel Body */}
            <div ref={scrollRef} className="flex-1 overflow-y-auto custom-scrollbar p-0">
                <AnimatePresence mode="wait">
                    {panelMode === 'CHAT' ? (
                        <motion.div
                            key="chat"
                            initial={{ opacity: 0 }}
                            animate={{ opacity: 1 }}
                            exit={{ opacity: 0 }}
                            className="p-4 space-y-4"
                        >
                            {/* Empty State / Intro */}
                            {logs.length === 0 && (
                                <div className="py-8 text-center space-y-3 opacity-60">
                                    <MessageSquare size={24} className="mx-auto text-gray-400" />
                                    <div className="space-y-1">
                                        <p className="text-[11px] font-bold text-black uppercase tracking-wider">Awaiting Instruction</p>
                                        <p className="text-[10px] text-gray-500">Describe a molecule or structural transformation.</p>
                                    </div>
                                </div>
                            )}

                            {/* Suggestions (Internalized) */}
                            {logs.length === 0 && (
                                <div className="space-y-2">
                                    {suggestions.map((s) => (
                                        <button
                                            key={s.label}
                                            onClick={() => setCommand(s.label)}
                                            className="w-full flex items-center gap-3 p-3 bg-white hover:bg-[#F9FAFB] border border-[#E5E7EB] rounded-lg text-left transition-all shadow-sm group"
                                        >
                                            <span className="text-[10px] font-mono text-blue-600/80 w-6 text-center group-hover:scale-110 transition-transform">{s.icon}</span>
                                            <span className="text-[11px] font-medium text-gray-600">{s.label}</span>
                                        </button>
                                    ))}
                                </div>
                            )}

                            {/* Brief Chat List */}
                            {logs.map((log) => (
                                <div key={log.id} className={`group border-l-2 transition-colors pl-4 py-2 ${log.source === 'system' ? 'border-red-200 bg-red-50/30' : 'border-[#E5E7EB] hover:border-blue-400'}`}>
                                    <div className={`text-[11px] font-medium leading-relaxed ${log.source === 'system' ? 'text-red-700' : 'text-black'}`}>
                                        {log.description}
                                    </div>
                                    <div className="flex items-center gap-2 mt-2">
                                        <span className={`text-[8px] font-black uppercase tracking-widest px-1.5 py-0.5 rounded-sm ${log.source === 'user' ? 'bg-blue-100 text-blue-700' : 'bg-gray-100 text-gray-500'}`}>
                                            {log.source === 'system' ? '⚠ System' : log.source}
                                        </span>
                                        <span className="text-[8px] text-gray-400 font-mono">
                                            {new Date(log.timestamp).toLocaleTimeString([], { hour: '2-digit', minute: '2-digit' })}
                                        </span>
                                    </div>
                                </div>
                            ))}
                        </motion.div>
                    ) : (
                        <motion.div
                            key="log"
                            initial={{ opacity: 0 }}
                            animate={{ opacity: 1 }}
                            exit={{ opacity: 0 }}
                            className="p-3 space-y-1"
                        >
                            {logs.length === 0 ? (
                                <div className="py-8 text-center opacity-40 space-y-2 font-mono">
                                    <p className="text-[10px] uppercase tracking-widest items-center flex justify-center gap-2 text-gray-400">
                                        No Runtime Data
                                    </p>
                                </div>
                            ) : (
                                logs.map((log) => (
                                    <div key={log.id} className="group flex items-start gap-3 p-2 border-b border-gray-100 hover:bg-gray-50 transition-colors">
                                        <div className="shrink-0 mt-1">
                                            <div className={`w-1 h-1 rounded-full ${log.source === 'system' ? 'bg-red-500' : 'bg-blue-500'}`} />
                                        </div>
                                        <div className="flex-1 min-w-0">
                                            <div className="flex justify-between items-center mb-0.5">
                                                <span className="text-[9px] font-mono text-gray-400 uppercase tabular-nums">
                                                    {new Date(log.timestamp).toLocaleTimeString([], { hour: '2-digit', minute: '2-digit', second: '2-digit' })}
                                                </span>
                                                <span className="text-[8px] font-bold text-gray-400 uppercase tracking-tighter">{log.mode}</span>
                                            </div>
                                            <p className={`text-[10px] font-mono break-words leading-tight ${log.source === 'system' ? 'text-red-600' : 'text-gray-600'}`}>
                                                {log.description}
                                            </p>
                                        </div>
                                    </div>
                                ))
                            )}
                        </motion.div>
                    )}
                </AnimatePresence>
            </div>

            {/* Panel Footer (Chat Input) */}
            <div className="p-4 bg-[#F9FAFB] border-t border-[#E5E7EB]">
                <form onSubmit={handleSubmit} className="relative group">
                    <input
                        type="text"
                        value={command}
                        onChange={(e) => setCommand(e.target.value)}
                        placeholder={isCommandRunning ? "Building..." : "Ask Gemini..."}
                        disabled={isCommandRunning}
                        className={`w-full bg-white border border-[#E5E7EB] ${isCommandRunning ? 'border-blue-400 ring-4 ring-blue-50' : 'hover:border-gray-300 focus:border-blue-500 focus:ring-1 focus:ring-blue-500/10'} 
                        rounded-lg px-4 py-3 text-[11px] outline-none transition-all placeholder:text-gray-400 pr-10 text-black shadow-sm`}
                    />
                    <button
                        type="submit"
                        disabled={!command.trim() || isCommandRunning}
                        className="absolute right-1 top-1/2 -translate-y-1/2 p-2 rounded-md hover:text-blue-600 text-gray-400 disabled:opacity-0 transition-all"
                    >
                        {isCommandRunning ? <Loader2 size={14} className="animate-spin text-blue-600" /> : <Send size={14} />}
                    </button>
                </form>
                <div className="mt-4 flex items-center justify-between text-[8px] font-black text-gray-400 uppercase tracking-[0.2em] px-1 transition-colors">
                    <div className="flex items-center gap-2">
                        <span className="w-1.5 h-1.5 rounded-full bg-green-500" />
                        <span>Ready</span>
                    </div>
                    <div className="flex flex-col items-end gap-0.5">
                        <span>L4-Model Node</span>
                        <span className="text-[6px] tracking-normal text-blue-500/60 lowercase italic font-medium">Chemistry engine: RDKit Deterministic</span>
                    </div>
                </div>
            </div>
        </div>
    );
}
