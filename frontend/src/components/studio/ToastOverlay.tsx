
import { useEffect, useState } from 'react';
import { useHistoryStore } from '../../store/historyStore';
import { motion, AnimatePresence } from 'framer-motion';
import { CheckCircle2 } from 'lucide-react';

export default function ToastOverlay() {
    const logs = useHistoryStore(state => state.logs);
    const [visibleToast, setVisibleToast] = useState<string | null>(null);
    const [lastLogId, setLastLogId] = useState<string | null>(null);

    useEffect(() => {
        if (logs.length > 0) {
            const latest = logs[logs.length - 1];
            if (latest.id !== lastLogId && latest.source !== 'system') {
                setLastLogId(latest.id);
                setVisibleToast(latest.description);
                const timer = setTimeout(() => setVisibleToast(null), 3000);
                return () => clearTimeout(timer);
            }
        }
    }, [logs, lastLogId]);

    return (
        <div className="absolute top-8 left-1/2 -translate-x-1/2 z-[100] pointer-events-none">
            <AnimatePresence>
                {visibleToast && (
                    <motion.div
                        initial={{ opacity: 0, y: -20, scale: 0.95 }}
                        animate={{ opacity: 1, y: 0, scale: 1 }}
                        exit={{ opacity: 0, scale: 0.95 }}
                        className="flex items-center gap-3 px-5 py-3 bg-black/90 backdrop-blur-xl border border-white/20 rounded-2xl shadow-2xl text-white"
                    >
                        <div className="p-1 bg-green-500/20 rounded-lg">
                            <CheckCircle2 size={16} className="text-green-500" />
                        </div>
                        <div className="flex flex-col">
                            <span className="text-[10px] uppercase font-black tracking-widest text-white/40 leading-none mb-1">State Updated</span>
                            <span className="text-xs font-bold leading-none">{visibleToast}</span>
                        </div>
                    </motion.div>
                )}
            </AnimatePresence>
        </div>
    );
}
