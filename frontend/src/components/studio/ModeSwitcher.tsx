
import { useStudioStore } from '../../store/studioStore';
import type { StudioMode } from '../../types/studio';

const MODES: StudioMode[] = ['design', 'optimize', 'simulate'];

export default function ModeSwitcher() {
    const { mode, setMode } = useStudioStore();

    return (
        <div className="bg-white p-1.5 rounded-xl border border-lightGrey shadow-sm flex items-center gap-1">
            {MODES.map((m) => (
                <button
                    key={m}
                    onClick={() => setMode(m)}
                    className={`px-6 py-2 rounded-lg text-sm font-bold transition-all ${mode === m
                            ? 'bg-black text-white shadow-md transform scale-105'
                            : 'text-midGrey hover:text-black hover:bg-offwhite'
                        }`}
                >
                    {m.charAt(0).toUpperCase() + m.slice(1)}
                </button>
            ))}
        </div>
    );
}
