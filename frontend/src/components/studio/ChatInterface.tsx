
import React, { useState, useRef, useEffect } from 'react';
import type { StudioMessage } from '../../types/studio';
import { StudioGateway } from '../../lib/ai/studioGateway';
import ChatMessage from './ChatMessage';
import { useStudioStore } from '../../store/studioStore';
import { moleculeToJSON } from '../../lib/engineAdapter';

export default function ChatInterface() {
    const { mode, molecule, messages, addMessage } = useStudioStore();
    const [inputValue, setInputValue] = useState('');
    const [isLoading, setIsLoading] = useState(false);
    const scrollRef = useRef<HTMLDivElement>(null);

    // Auto-scroll to bottom
    useEffect(() => {
        if (scrollRef.current) {
            scrollRef.current.scrollTop = scrollRef.current.scrollHeight;
        }
    }, [messages]);

    const handleSubmit = async (e?: React.FormEvent) => {
        e?.preventDefault();
        if (!inputValue.trim() || isLoading) return;

        const userMessage: StudioMessage = {
            id: Date.now().toString(),
            role: 'user',
            content: inputValue,
            timestamp: new Date(),
            metadata: {
                mode,
            }
        };

        addMessage(userMessage);
        setInputValue('');
        setIsLoading(true);

        try {
            // Serialize molecule context.
            // If molecule is null (shouldn't be in main flow), we pass empty struct or null
            const moleculeContext = molecule ? moleculeToJSON(molecule) : undefined;

            const response = await StudioGateway.chat({
                prompt: userMessage.content,
                mode,
                moleculeContext,
                history: messages,
            });

            const mentorMessage: StudioMessage = {
                id: (Date.now() + 1).toString(),
                role: 'mentor',
                content: response.response,
                timestamp: new Date(),
                metadata: {
                    mode,
                }
            };

            addMessage(mentorMessage);
        } catch (error) {
            console.error("ChatInterface Error:", error);
            const errorMessage: StudioMessage = {
                id: (Date.now() + 1).toString(),
                role: 'mentor',
                content: "System Error: Unable to process command.",
                timestamp: new Date(),
                metadata: { mode }
            };
            addMessage(errorMessage);
        } finally {
            setIsLoading(false);
        }
    };

    return (
        <div className="flex flex-col h-full bg-white">
            {/* Chat Messages */}
            <div
                ref={scrollRef}
                className="flex-1 overflow-y-auto px-4 py-4 scrollbar-hide space-y-4"
            >
                {messages.length === 0 && (
                    <div className="flex flex-col items-center justify-center h-full opacity-40 text-center px-8">
                        <div className="w-12 h-12 border border-black/10 rounded-full flex items-center justify-center mb-4 text-black">
                            {mode === 'design' && <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2"><path d="M12 20h9" /><path d="M16.5 3.5a2.121 2.121 0 0 1 3 3L7 19l-4 1 1-4L16.5 3.5z" /></svg>}
                            {mode === 'optimize' && <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2"><circle cx="12" cy="12" r="10" /><line x1="12" y1="16" x2="12" y2="12" /><line x1="12" y1="8" x2="12.01" y2="8" /></svg>}
                            {mode === 'simulate' && <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2"><path d="M2 12h20" /><path d="M2 12l5-5" /><path d="M2 12l5 5" /></svg>}
                        </div>
                        <p className="text-sm font-medium text-black">
                            {mode === 'design' ? 'Ready to design.' :
                                mode === 'optimize' ? 'Ready to optimize.' :
                                    'Ready to simulate.'}
                        </p>
                        <p className="text-xs mt-1 text-midGrey">Type a command to begin.</p>
                    </div>
                )}

                {messages.map((msg) => (
                    <ChatMessage key={msg.id} message={msg} />
                ))}

                {isLoading && (
                    <div className="flex justify-start mb-4 animate-pulse">
                        <div className="flex gap-1 bg-offwhite px-3 py-2 rounded-lg">
                            <div className="w-1.5 h-1.5 bg-midGrey rounded-full animate-bounce" style={{ animationDelay: '0ms' }} />
                            <div className="w-1.5 h-1.5 bg-midGrey rounded-full animate-bounce" style={{ animationDelay: '150ms' }} />
                            <div className="w-1.5 h-1.5 bg-midGrey rounded-full animate-bounce" style={{ animationDelay: '300ms' }} />
                        </div>
                    </div>
                )}
            </div>

            {/* Input Area */}
            <form onSubmit={handleSubmit} className="border-t border-lightGrey p-4 bg-white relative">
                <input
                    type="text"
                    value={inputValue}
                    onChange={(e) => setInputValue(e.target.value)}
                    placeholder={
                        mode === 'design' ? "Add a methyl group..." :
                            mode === 'optimize' ? "Improve solubility..." :
                                "Run MD simulation..."
                    }
                    disabled={isLoading}
                    className="w-full bg-offwhite border border-lightGrey rounded-xl px-4 py-3.5 pr-12 text-sm text-black placeholder:text-midGrey focus:outline-none focus:ring-2 focus:ring-black/5 focus:border-black/20 transition-all disabled:opacity-50"
                />
                <button
                    type="submit"
                    disabled={isLoading || !inputValue.trim()}
                    className="absolute right-6 top-1/2 -translate-y-1/2 text-black hover:text-blue-600 transition-colors disabled:text-lightGrey"
                >
                    <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
                        <line x1="22" y1="2" x2="11" y2="13"></line>
                        <polygon points="22 2 15 22 11 13 2 9 22 2"></polygon>
                    </svg>
                </button>
            </form>
        </div>
    );
}
