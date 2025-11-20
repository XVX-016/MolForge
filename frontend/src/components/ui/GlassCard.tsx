import React from 'react';

interface GlassCardProps {
  children: React.ReactNode;
  className?: string;
}

export default function GlassCard({ children, className = '' }: GlassCardProps) {
  return (
    <div
      className={
        'backdrop-blur-xl bg-white/10 border border-white/20 rounded-2xl shadow-lg p-4 ' +
        className
      }
    >
      {children}
    </div>
  );
}

