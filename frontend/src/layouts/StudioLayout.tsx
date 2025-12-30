import { Link, Outlet } from 'react-router-dom';
import { Box, User } from 'lucide-react';

export default function StudioLayout() {
    return (
        <div className="flex flex-col h-screen bg-[#0B0B0D] text-[#E5E7EB] font-sans selection:bg-blue-500/30 overflow-hidden">
            {/* Discrete Studio Header */}
            <header className="h-14 shrink-0 px-6 border-b border-[#1F2937] bg-[#0F1013] flex items-center justify-between z-50">
                <div className="flex items-center gap-6">
                    <Link to="/" className="flex items-center gap-2 group">
                        <div className="w-8 h-8 bg-blue-600 rounded-lg flex items-center justify-center shadow-[0_0_15px_rgba(37,99,235,0.4)] group-hover:scale-105 transition-transform">
                            <Box size={18} className="text-white" />
                        </div>
                        <span className="text-sm font-black uppercase tracking-tighter text-white">
                            MolForge <span className="text-blue-500">Studio</span>
                        </span>
                    </Link>

                    <nav className="flex items-center gap-1 ml-4">
                        <Link to="/" className="px-3 py-1 text-[11px] font-bold text-gray-500 hover:text-white transition-colors uppercase tracking-widest">Home</Link>
                        <Link to="/library" className="px-3 py-1 text-[11px] font-bold text-gray-500 hover:text-white transition-colors uppercase tracking-widest">Library</Link>
                        <Link to="/docs" className="px-3 py-1 text-[11px] font-bold text-gray-500 hover:text-white transition-colors uppercase tracking-widest">Docs</Link>
                    </nav>
                </div>

                <div className="flex items-center gap-4">
                    <div className="flex items-center gap-2 px-3 py-1 bg-white/5 border border-white/10 rounded-full">
                        <div className="w-1.5 h-1.5 rounded-full bg-green-500 animate-pulse" />
                        <span className="text-[9px] font-black uppercase tracking-widest text-gray-400">Environment: Stable</span>
                    </div>

                    <button className="w-8 h-8 rounded-full bg-white/5 border border-white/10 flex items-center justify-center hover:bg-white/10 transition-colors">
                        <User size={16} className="text-gray-400" />
                    </button>
                </div>
            </header>

            {/* Studio Content */}
            <main className="flex-1 overflow-hidden">
                <Outlet />
            </main>
        </div>
    );
}
