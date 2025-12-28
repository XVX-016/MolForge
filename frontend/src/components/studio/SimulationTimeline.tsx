
import React, { useState, useEffect, useRef } from 'react';
import { Play, Pause, RotateCcw } from 'lucide-react';
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
            // Throttling to approx 30fps for visibility
            const delta = time - lastTimeRef.current;
            if (delta > 33) { // ~30ms
                setCurrentStep(prev => {
                    const next = prev + 1;
                    if (next >= trajectory.length) return 0; // Loop
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

    // Notify parent of frame change
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
        <div className="absolute inset-x-0 bottom-0 bg-white border-t border-lightGrey p-4 flex flex-col gap-2 z-30">
            {/* Header / Info */}
            <div className="flex justify-between items-end">
                <div className="flex flex-col">
                    <span className="text-[10px] uppercase font-bold text-midGrey tracking-wider">Simulation Time</span>
                    <span className="text-xl font-mono font-bold text-black">
                        {currentFrame.time.toFixed(1)} <span className="text-sm font-normal text-midGrey">ps</span>
                    </span>
                </div>
                <div className="flex flex-col items-end">
                    <span className="text-[10px] uppercase font-bold text-midGrey tracking-wider">Potential Energy</span>
                    <span className="text-sm font-mono font-bold text-black">
                        {currentFrame.energy.toFixed(2)} <span className="text-xs font-normal text-midGrey">kcal/mol</span>
                    </span>
                </div>
            </div>

            {/* Controls Row */}
            <div className="flex items-center gap-4">

                <div className="flex items-center gap-2">
                    <button
                        onClick={() => setCurrentStep(0)}
                        className="p-2 rounded-full hover:bg-gray-100 text-darkGrey transition-colors"
                        title="Reset"
                    >
                        <RotateCcw size={16} />
                    </button>
                    <button
                        onClick={() => setIsPlaying(!isPlaying)}
                        className="p-3 rounded-full bg-black text-white hover:bg-gray-800 transition-colors shadow-md"
                        title={isPlaying ? "Pause" : "Play"}
                    >
                        {isPlaying ? <Pause size={20} fill="currentColor" /> : <Play size={20} fill="currentColor" className="ml-0.5" />}
                    </button>
                </div>

                {/* Scrubber */}
                <div className="flex-1 relative h-8 flex items-center">
                    <input
                        type="range"
                        min={0}
                        max={trajectory.length - 1}
                        value={currentStep}
                        onChange={handleScrub}
                        className="w-full h-2 bg-gray-200 rounded-lg appearance-none cursor-pointer accent-black"
                    />
                </div>

                <span className="font-mono text-xs text-midGrey w-12 text-right">
                    {currentStep}/{trajectory.length - 1}
                </span>

            </div>
        </div>
    );
}
