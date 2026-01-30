import StudioTopBar from '../components/studio/StudioTopBar';
import MetricsBar from '../components/studio/MetricsBar';
import WorkflowTimeline from '../components/studio/WorkflowTimeline';
import MultimodalWorkspace from '../components/studio/MultimodalWorkspace';
import GeminiDesk from '../components/studio/GeminiDesk';
import { useStudioV2Store } from '../store/useStudioV2Store';

export default function StudioPage() {
  const {
    currentExperiment,
    nodes,
    selectedNodeId,
    compareNodeId,
    activeInsight,
    loading,
    addNode,
    runNode,
    selectNode,
    setCompareNode
  } = useStudioV2Store();

  // REMOVED: Auto-init effect.
  // The Studio now starts in an atomic "Empty" state, waiting for user action.

  const activeNode = nodes.find(n => n.id === selectedNodeId);

  return (
    <div className="flex flex-col h-screen overflow-hidden bg-[#F3F4F6]">
      <StudioTopBar />
      <MetricsBar />

      <div className="flex-1 flex overflow-hidden">
        {/* Left: Workflow Timeline */}
        <aside className="w-[320px] shrink-0 border-r border-[#E5E7EB] bg-white flex flex-col">
          <WorkflowTimeline
            nodes={nodes}
            selectedNodeId={selectedNodeId || undefined}
            compareNodeId={compareNodeId || undefined}
            onNodeSelect={selectNode}
            onNodeCompare={setCompareNode}
            onRunNode={runNode}
            onAddNode={addNode}
            isEmpty={!currentExperiment} // New prop
          />
        </aside>

        {/* Center: Multimodal Workspace */}
        <main className="flex-1 relative bg-[#F9FAFB] flex flex-col">
          <MultimodalWorkspace
            activeNode={activeNode}
            compareNode={nodes.find(n => n.id === compareNodeId)}
            onCloseComparison={() => setCompareNode(null)}
            _moleculeData={null}
            isEmpty={!currentExperiment} // New prop
          />
        </main>

        {/* Right: Gemini Reasoning Desk */}
        <aside className="w-[380px] shrink-0 border-l border-[#E5E7EB] bg-white flex flex-col overflow-y-auto">
          <GeminiDesk
            insight={activeInsight || undefined}
            loading={loading}
            isEmpty={!currentExperiment} // New prop
          />
        </aside>
      </div>
    </div>
  );
}
