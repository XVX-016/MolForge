import Navbar from "../Navbar";
import TopActionBar from "./layout/TopActionBar";
import LeftToolDock from "./layout/LeftToolDock";
import RightInspector from "./layout/RightInspector";
import BottomTemplateBar from "./layout/BottomTemplateBar";
import LabCanvas from "./LabCanvas";

export default function LabV2Page() {
    return (
        <div className="w-full h-screen flex flex-col overflow-hidden bg-gray-50 relative">
            {/* 0. Global Navbar (Likely staying, as Zone A starts at top-20 ~ 80px) */}
            <div className="z-[60] relative bg-white border-b border-gray-100">
                <Navbar />
            </div>

            {/* Main Content Area */}
            <div className="flex-1 relative w-full h-full">

                {/* Zone C: Canvas (Sacred) - z-0 */}
                <div id="lab-canvas" className="absolute inset-0 w-full h-full z-0 bg-[#F9FAFB]">
                    <LabCanvas />
                </div>

                {/* Zone A: Top Floating Action Bar - z-50 */}
                <TopActionBar />

                {/* Zone B: Left Tool Dock - z-50 */}
                <LeftToolDock />

                {/* Zone D: Right Inspector - z-50 */}
                <RightInspector />

                {/* Zone E: Bottom Template Bar - z-50 */}
                <BottomTemplateBar />

                {/* Event Blocker for background to let events pass to canvas where no UI exists */}
                {/* Actually, with fixed positioning of UI elements, we don't need a blocker if they don't cover everything. 
                    The canvas is at z-0, interactive. The UI is z-50, interactive. 
                    Clicks on 'empty' areas of UI layer (if it was a full div) would block. 
                    But here we use individual components with 'fixed' position. 
                    So clicks pass through to canvas naturally. */}
            </div>
        </div>
    );
}
