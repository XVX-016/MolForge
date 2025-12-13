# Lab Layout & Grid Refactor

## 🏗️ Phase 1: Layout Skeleton Implemented

As per the authoritative design goals, the **Lab Page** has been completely restructured into 5 distinct zones:

### **Zone A: Top Floating Action Bar**
*   **Location**: Fixed, Top-Center (`top-20`).
*   **Design**: Floating pill, icon-only, glassmorphism (`backdrop-blur`).
*   **Function**: Global actions (Select, Add, Undo/Redo, Save).

### **Zone B: Left Tool Dock**
*   **Location**: Fixed, Left (`left-4`, `top-28`).
*   **Design**: Collapsible sidebar. Shows icons by default, expands to show details on hover.
*   **Function**: Atom picker (color-coded), Bond types.

### **Zone C: The Sacred Canvas**
*   **Location**: Fullscreen underlying layer (`z-0`).
*   **Change**: **REMOVED VERTICAL GRIDLINES**. Now uses a custom horizontal-only grid system to keep the view clean and maximize uninterrupted space.

### **Zone D: Right Inspector**
*   **Location**: Fixed, Right (`right-4`).
*   **Design**: Sliding panel. Hidden by default (`translate-x-full`).
*   **Function**: Context-aware properties (placeholders for now).

### **Zone E: Bottom Template Bar**
*   **Location**: Fixed, Bottom-Center (`bottom-6`).
*   **Design**: Floating bar with drag-and-drop templates.

## 🎨 Styling
*   Implemented using Tailwind CSS specs provided.
*   Ensured UI elements float *above* the canvas (`z-50`) without blocking interaction on empty space.

## 🔜 Next Steps (Phase 2 & 3)
*   Wire up the **Right Inspector** to slide in when an atom is selected.
*   Connect the **Bottom Template Bar** to actual drag-and-drop logic.
*   Refine the interactions (hover states, animations).
