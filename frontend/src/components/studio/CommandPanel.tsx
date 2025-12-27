
import React, { useState, useRef, useEffect } from 'react';
import { useStudioStore } from '../../store/studioStore';
import { StudioGateway } from '../../lib/ai/studioGateway';
import type { StudioMessage } from '../../types/studio';
import { moleculeToJSON } from '../../lib/engineAdapter';
import ChatMessage from './ChatMessage';

export default function CommandPanel() {
    const { mode, molecule, messages, addMessage, setMolecule } = useStudioStore();
    const [inputValue, setInputValue] = useState('');
    const [isLoading, setIsLoading] = useState(false);
    const scrollRef = useRef<HTMLDivElement>(null);

    // Auto-scroll to bottom
    useEffect(() => {
        if (scrollRef.current) {
            scrollRef.current.scrollTop = scrollRef.current.scrollHeight;
        }
    }, [messages]);

    const handleCommand = async (e?: React.FormEvent) => {
        e?.preventDefault();
        if (!inputValue.trim() || isLoading) return;

        const commandText = inputValue;

        // 1. Add User Command to History
        const userMessage: StudioMessage = {
            id: Date.now().toString(),
            role: 'user',
            content: commandText,
            timestamp: new Date(),
            metadata: { mode }
        };
        addMessage(userMessage);
        setInputValue('');
        setIsLoading(true);

        try {
            // 2. Prepare Context
            const moleculeContext = molecule ? moleculeToJSON(molecule) : undefined;

            // 3. Execute Command via Gateway
            // Note: In a real implementation this would change to executeCommand signature
            const response = await StudioGateway.chat({
                prompt: commandText,
                mode,
                moleculeContext,
                history: messages,
            });

            // 4. Handle System Response
            const systemMessage: StudioMessage = {
                id: (Date.now() + 1).toString(),
                role: 'mentor', // Rename to 'system' in types later
                content: response.response,
                timestamp: new Date(),
                metadata: { mode }
            };
            addMessage(systemMessage);

            // TODO: In future, if response includes a patch, apply it:
            // if (response.moleculePatch) setMolecule(response.moleculePatch);

        } catch (error) {
            console.error("Command Error:", error);
            addMessage({
                id: Date.now().toString(),
                role: 'mentor',
                content: "System Error: Command failed execution.",
                timestamp: new Date(),
                metadata: { mode }
            });
        } finally {
            setIsLoading(false);
        }
    };

    return (
        <div className="flex flex-col h-full bg-white">
            {/* Header / Status Bar */}
            <div className="p-4 border-b border-lightGrey bg-offwhite/50 backdrop-blur-sm shrink-0">
                <div className="flex items-center justify-between">
                    <span className="text-xs font-bold uppercase tracking-wider text-midGrey">Command History</span>
                    <span className={`px-2 py-0.5 rounded text-[10px] font-bold uppercase ${mode === 'design' ? 'bg-blue-100 text-blue-700' :
                            mode === 'optimize' ? 'bg-purple-100 text-purple-700' :
                                'bg-orange-100 text-orange-700'
                        }`}>
                        {mode}
                    </span>
                </div>
            </div>

            {/* History Output */}
            <div
                ref={scrollRef}
                className="flex-1 overflow-y-auto px-4 py-4 scrollbar-hide space-y-4 font-mono text-sm"
            >
                {messages.length === 0 && (
                    <div className="flex flex-col items-center justify-center h-full opacity-40 text-center px-8 text-darkGrey">
                        <p className="font-semibold mb-2">
                            {mode === 'design' ? 'Design Mode Active' :
                                mode === 'optimize' ? 'Optimization Mode Active' :
                                    'Simulation Mode Active'}
                        </p>
                        <p className="text-xs">
                            {mode === 'design' && "Describe a structure or modification..."}
                            {mode === 'optimize' && "Describe a property to improve..."}
                            {mode === 'simulate' && "Describe conditions or reactions..."}
                        </p>
                    </div>
                )}

                {messages.map((msg) => (
                    <ChatMessage key={msg.id} message={msg} />
                ))}

                {isLoading && (
                    <div className="flex justify-start mb-4">
                        <span className="text-xs text-midGrey animate-pulse">Processing command...</span>
                    </div>
                )}
            </div>

            {/* Command Input */}
            <form onSubmit={handleCommand} className="border-t border-lightGrey p-4 bg-white relative shrink-0">
                <div className="flex items-center gap-2 bg-offwhite border border-lightGrey rounded-xl px-4 py-3 focus-within:ring-2 focus-within:ring-black/5 focus-within:border-black/20 transition-all">
                    <span className="text-midGrey font-mono select-none">{'>'}</span>
                    <input
                        type="text"
                        value={inputValue}
                        onChange={(e) => setInputValue(e.target.value)}
                        placeholder="Type a command..."
                        disabled={isLoading}
                        className="flex-1 bg-transparent border-none p-0 text-sm text-black placeholder:text-midGrey focus:outline-none font-mono"
                        autoFocus
                    />
                </div>
            </form>
        </div>
    );
}
