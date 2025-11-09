# Frontend UI System Documentation

## Overview

The frontend UI system provides a complete molecular builder interface with tool-based editing, undo/redo, and real-time visualization.

## Tool System Architecture

### Tool Modes

**Location:** `frontend/src/store/moleculeStore.ts`

**Tool Types:**
- `select` - Default mode for selecting and dragging atoms
- `add-atom` - Place new atoms by clicking empty space
- `bond` - Create bonds between atoms (single/double/triple)
- `delete` - Delete atoms or bonds

**State Management:**
- `tool: Tool` - Current active tool
- `atomToAdd: Element | null` - Element selected for placement
- `currentBondOrder: number` - Bond order (1, 2, or 3)

**Actions:**
- `setTool(tool)` - Switch tool mode
- `resetTool()` - Reset to select mode
- `setAtomToAdd(element)` - Set element for atom placement
- `setBondOrder(order)` - Set bond order

### ToolPanel Component

**Location:** `frontend/src/components/ToolPanel.tsx`

**Features:**
- Floating left toolbar
- Four tool buttons with icons
- Bond order selector (when bond tool active)
- Undo/Redo buttons
- Keyboard shortcuts (S, A, B, DEL)
- Visual highlighting of active tool

**Keyboard Shortcuts:**
- `S` - Select tool
- `A` - Add atom tool
- `B` - Bond tool
- `DEL` / `Backspace` - Delete tool
- `Ctrl+Z` - Undo
- `Ctrl+Y` / `Ctrl+Shift+Z` - Redo

### AtomPalette Component

**Location:** `frontend/src/components/AtomPalette.tsx`

**Features:**
- Grid of 10 elements (H, C, N, O, F, S, P, Cl, Br, I)
- Color-coded by element
- Auto-shows when add-atom tool is active
- Clicking element sets `atomToAdd` and switches to add-atom tool

## Atom Placement

### Add Atom Tool Flow

1. User clicks "Add Atom" button → `setTool('add-atom')`
2. AtomPalette appears
3. User selects element → `setAtomToAdd(element)`
4. User clicks empty space in canvas
5. `CanvasClickHandler` converts click to world coordinates
6. `addAtom(element, position)` called
7. Atom added to MoleculeGraph
8. Geometry optimized (5 iterations)
9. Predictions re-run

### Implementation

**CanvasClickHandler** (`frontend/src/components/r3f/CanvasClickHandler.tsx`):
- Listens for canvas clicks when `tool === 'add-atom'`
- Converts screen coordinates to 3D world coordinates
- Calls `engineAdapter.addAtom()`

**engineAdapter.addAtom()**:
- Creates new molecule if none exists
- Adds atom to MoleculeGraph
- Pushes state to history
- Triggers geometry optimization
- Updates store and triggers predictions

## Bond Tools

### Bond Creation Flow

1. User clicks "Bond" button → `setTool('bond')`
2. User selects bond order (Single/Double/Triple)
3. User clicks first atom → stored in `firstSelectedRef`
4. User clicks second atom → `addBond(a1, a2, order)` called
5. Bond created with specified order
6. Geometry optimized (10 iterations)
7. Predictions re-run

### Bond Order

**Supported Orders:**
- **Single (1)** - Default, standard covalent bond
- **Double (2)** - Double bond (e.g., C=O)
- **Triple (3)** - Triple bond (e.g., N≡N)

**Visual Representation:**
- Single: radius 0.14
- Double: radius 0.18
- Triple: radius 0.22

### BondTool Hook

**Location:** `frontend/src/components/r3f/BondTool.ts`

**Logic:**
- Only active when `tool === 'bond'`
- Tracks first selected atom
- Creates bond when second atom selected
- Prevents duplicate bonds
- Uses `currentBondOrder` from store

## Delete Tool

### Delete Flow

1. User clicks "Delete" button → `setTool('delete')`
2. Atoms/bonds show red outline on hover
3. User clicks atom or bond
4. `removeAtom(id)` or `removeBond(id)` called
5. Item removed from MoleculeGraph
6. Connected bonds removed (if atom deleted)
7. Selection cleared if deleted item was selected
8. Predictions re-run

### Visual Feedback

- **Hover:** Red outline glow
- **Selected:** Red outline + scale
- **Cursor:** `not-allowed` when hovering deletable items

### Implementation

**AtomMesh:**
- Checks `tool === 'delete'` on click
- Calls `removeAtom(id)` dynamically

**BondMesh:**
- Checks `tool === 'delete'` on click
- Calls `removeBond(id)` dynamically

**engineAdapter:**
- `removeAtom()` - Removes atom and all connected bonds
- `removeBond()` - Removes single bond
- Both push state to history

## Undo/Redo System

### History Store

**Location:** `frontend/src/store/historyStore.ts`

**Architecture:**
- Uses `UndoStack` from engine package
- Stores MoleculeGraph snapshots
- Tracks `canUndo` and `canRedo` states

**Actions:**
- `pushState()` - Save current molecule state
- `undo()` - Restore previous state
- `redo()` - Restore next state
- `clearHistory()` - Clear all history

### Integration

**Automatic State Saving:**
- `addAtom()` → `pushState()`
- `removeAtom()` → `pushState()`
- `removeBond()` → `pushState()`
- `addBond()` → `pushState()`

**Manual Undo/Redo:**
- ToolPanel buttons
- Keyboard shortcuts (Ctrl+Z, Ctrl+Y)

### State Serialization

- MoleculeGraph cloned before pushing
- Uses `MoleculeGraph.clone()` for deep copy
- Restores complete molecule structure

## Selection & Hover Enhancements

### Visual Feedback by Tool

**Select Tool:**
- Hover: Blue outline
- Selected: Blue outline + 1.15x scale

**Bond Tool:**
- Hover: Blue outline
- Selected: Blue outline + 1.15x scale
- Indicates atoms available for bonding

**Delete Tool:**
- Hover: Red outline
- Selected: Red outline + 1.15x scale
- Warning that item will be deleted

**Add Atom Tool:**
- Cursor: `crosshair`
- Ghost preview (TODO: implement)

### Implementation

**AtomMesh:**
- Reads `tool` from store
- Adjusts `outlineColor` based on tool
- Red for delete, blue for others

**BondMesh:**
- Similar tool-based coloring
- Red outline for delete tool

## Property Badges

### Atom Properties

**Displayed:**
- Element symbol
- 3D position (x, y, z)
- Atom ID

**Location:** Bottom-left of MoleculeViewer

**Updates:** Automatically when atom selected

### Bond Properties

**Displayed:**
- Bond order (Single/Double/Triple)
- Bond length (Å)
- Connected atoms (element symbols)

**Location:** Bottom-left of MoleculeViewer

**Updates:** Automatically when bond selected

## Rendering Pipeline

1. **MoleculeGraph** → `moleculeToRenderable()` → Renderable format
2. **AtomMesh** components render with tool-aware visuals
3. **BondMesh** components render with order-based radius
4. **ToolPanel** provides tool selection
5. **AtomPalette** shows when add-atom active
6. **PropertyBadge** shows selection details
7. **CanvasClickHandler** handles empty space clicks

## Event Flow

```
User Action → ToolPanel/Tool → Store Update
                ↓
         MoleculeViewer Reacts
                ↓
         AtomMesh/BondMesh Update Visuals
                ↓
         Click Handler (if applicable)
                ↓
         engineAdapter Function
                ↓
         pushState() → History
                ↓
         MoleculeGraph Update
                ↓
         Store.setMolecule()
                ↓
         Re-render + fetchPredictions()
```

## Keyboard Shortcuts

| Key | Action |
|-----|--------|
| `S` | Select tool |
| `A` | Add atom tool |
| `B` | Bond tool |
| `DEL` / `Backspace` | Delete tool |
| `Ctrl+Z` | Undo |
| `Ctrl+Y` | Redo |
| `Escape` | Cancel bond creation |

## Testing

**Unit Tests:**
- `ToolPanel.test.tsx` - Tool selection
- `historyStore.test.ts` - Undo/redo functionality
- `engineAdapter.addAtom.test.ts` - Atom operations

## Future Enhancements

- [ ] Ghost atom preview in add-atom mode
- [ ] Temporary line preview in bond mode
- [ ] Multi-select for batch operations
- [ ] Copy/paste molecules
- [ ] Tool presets (common molecule templates)
- [ ] Measurement tools (distance, angle)
- [ ] Valence validation warnings

