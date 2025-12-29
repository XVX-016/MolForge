
import React, { useState, useEffect, useRef } from 'react';
import { Play, Pause, RotateCcw, Activity } from 'lucide-react';
import type { Trajectory, SimulationFrame } from '../../lib/simulation';

interface SimulationTimelineProps {
    trajectory: Trajectory;
    onFrameChange: (frame: SimulationFrame) => void;
}

export default function SimulationTimeline({ trajectory, onFrameChange }: SimulationTimelineProps) {
    const [isPlaying, setIsPlaying] = useState(false);
    const [currentStep, setCurrentStep] = useState(0);
    const requestRef = useRef<number>();
    const lastTimeRef = useRef<number>();

    // Derived current frame
    const currentFrame = trajectory[currentStep];

    const animate = (time: number) => {
        if (lastTimeRef.current !== undefined) {
            const delta = time - lastTimeRef.current;
            if (delta > 33) { // ~30fps
                setCurrentStep(prev => {
                    const next = prev + 1;
                    if (next >= trajectory.length) return 0;
                    return next;
                });
                lastTimeRef.current = time;
            }
        } else {
            lastTimeRef.current = time;
        }
        requestRef.current = requestAnimationFrame(animate);
    };

    useEffect(() => {
        if (isPlaying) {
            requestRef.current = requestAnimationFrame(animate);
        } else {
            if (requestRef.current) cancelAnimationFrame(requestRef.current);
            lastTimeRef.current = undefined;
        }
        return () => {
            if (requestRef.current) cancelAnimationFrame(requestRef.current);
        };
    }, [isPlaying, trajectory.length]);

    useEffect(() => {
        if (currentFrame) {
            onFrameChange(currentFrame);
        }
    }, [currentStep, currentFrame, onFrameChange]);

    const handleScrub = (e: React.ChangeEvent<HTMLInputElement>) => {
        setIsPlaying(false);
        setCurrentStep(Number(e.target.value));
    };

    return (
        <div className="absolute inset-x-6 bottom-6 flex flex-col gap-4 z-30 animate-in fade-in slide-in-from-bottom-4 duration-500">
            {/* Top Stats Bar */}
            <div className="flex justify-between items-center px-6 py-3 bg-black/90 backdrop-blur-xl border border-white/10 rounded-2xl shadow-2xl">
                <div className="flex items-center gap-6">
                    <div className="flex flex-col">
                        <span className="text-[9px] uppercase font-black text-white/40 tracking-[0.2em] mb-0.5">Dynamics Clock</span>
                        <span className="text-lg font-mono font-bold text-white leading-none">
                            {currentFrame.time.toFixed(1)} <span className="text-xs font-normal text-white/60">ps</span>
                        </span>
                    </div>
                    <div className="w-px h-8 bg-white/10" />
                    <div className="flex flex-col">
                        <span className="text-[9px] uppercase font-black text-white/40 tracking-[0.2em] mb-0.5">System Energy</span>
                        <span className="text-sm font-mono font-bold text-orange-400 leading-none">
                            {currentFrame.energy.toFixed(2)} <span className="text-[10px] font-normal text-white/60">kcal/mol</span>
                        </span>
                    </div>
                </div>

                <div className="flex items-center gap-2 px-3 py-1 bg-white/5 rounded-full border border-white/5">
                    <Activity size={12} className="text-orange-500" />
                    <span className="text-[10px] font-bold text-white/80 uppercase tracking-widest">Active Simulation</span>
                </div>
            </div>

            {/* Main Interactive Timeline */}
            <div className="bg-white/95 backdrop-blur-md border border-lightGrey/50 rounded-2xl p-5 shadow-2xl flex flex-col gap-5">
                <div className="flex items-center justify-between">
                    <div className="flex items-center gap-5">
                        <button
                            onClick={() => setIsPlaying(!isPlaying)}
                            className={`w-14 h-14 flex items-center justify-center rounded-2xl transition-all shadow-xl active:scale-95 ${isPlaying
                                ? 'bg-orange-500 text-white shadow-orange-500/20'
                                : 'bg-black text-white hover:bg-gray-800'
                                }`}
                        >
                            {isPlaying ? <Pause size={24} fill="currentColor" /> : <Play size={24} fill="currentColor" className="ml-1" />}
                        </button>

                        <div className="flex flex-col">
                            <h3 className="text-sm font-black text-black leading-tight">Vibrational Replay</h3>
                            <p className="text-[11px] font-medium text-midGrey">Deterministic trajectory analysis</p>
                        </div>
                    </div>

                    <div className="flex items-center gap-2">
                        <button
                            onClick={() => setCurrentStep(0)}
                            className="p-3 rounded-xl hover:bg-gray-100 text-darkGrey transition-all border border-transparent hover:border-lightGrey/50"
                            title="Reset to 0ps"
                        >
                            <RotateCcw size={18} />
                        </button>
                    </div>
                </div>

                {/* Scrubber + Markers */}
                <div className="flex flex-col gap-2">
                    <div className="relative h-6 flex items-center group">
                        <div className="absolute inset-0 h-1.5 bg-gray-100 rounded-full my-auto" />
                        <input
                            type="range"
                            min={0}
                            max={trajectory.length - 1}
                            value={currentStep}
                            onChange={handleScrub}
                            className="absolute inset-0 w-full h-1.5 bg-transparent appearance-none cursor-pointer accent-orange-500 z-10"
                        />
                        {/* Interactive progress bar */}
                        <div
                            className="absolute left-0 h-1.5 bg-orange-500 rounded-full my-auto pointer-events-none transition-all duration-100"
                            style={{ width: `${(currentStep / (trajectory.length - 1)) * 100}%` }}
                        />
                    </div>

                    <div className="flex justify-between text-[10px] font-black text-gray-400 uppercase tracking-tighter">
                        <span>0.0 ps</span>
                        <span className="text-black bg-gray-100 px-2 rounded-full">Step {currentStep} / {trajectory.length - 1}</span>
                        <span>{(trajectory.length * 0.5).toFixed(1)} ps</span>
                    </div>
                </div>
            </div>
        </div>
    );
}
