import { useEffect } from 'react';
import { useStudioStore } from '../store/studioStore';
import StudioTopBar from '../components/studio/StudioTopBar';
import MetricsBar from '../components/studio/MetricsBar';
import IntentPanel from '../components/studio/StudioControlPanel';
import StudioMainCanvas from '../components/studio/StudioMainCanvas';
import AuditPanel from '../components/studio/StudioPropertiesPanel';

export default function StudioPage() {
  const { loadDashboard } = useStudioStore();

  useEffect(() => {
    // Initial load - in a real app, we'd get the ID from the URL/Library
    // Hardcoding a placeholder for now to trigger the READY state
    loadDashboard("00000000-0000-0000-0000-000000000000");
  }, [loadDashboard]);

  return (
    <div className="flex flex-col h-screen overflow-hidden bg-[#F3F4F6]">
      <StudioTopBar />
      <MetricsBar />

      <div className="flex-1 flex overflow-hidden">
        {/* Left: Intent & Control */}
        <aside className="w-[380px] shrink-0 border-r border-[#E5E7EB] bg-white flex flex-col">
          <IntentPanel />
        </aside>

        {/* Center: Structural Truth */}
        <main className="flex-1 relative bg-[#F9FAFB] flex flex-col">
          <StudioMainCanvas />
        </main>

        {/* Right: Clinical Reality */}
        <aside className="w-[400px] shrink-0 border-l border-[#E5E7EB] bg-white flex flex-col overflow-y-auto">
          <AuditPanel />
        </aside>
      </div>
    </div>
  );
}
