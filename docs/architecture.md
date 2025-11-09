# BioSynth AI Architecture

## Overview

BioSynth AI is a full-stack molecular design system with three main components:

1. **Frontend** - React + TypeScript + R3F for 3D visualization
2. **Engine** - Pure TypeScript molecular graph engine
3. **Backend** - FastAPI + PyTorch for ML predictions

## Architecture Diagram

```
┌─────────────────────────────────────────────────────────────┐
│                         Frontend                             │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐       │
│  │   Dashboard  │  │ MoleculeView│  │   Store      │       │
│  │              │  │    er        │  │  (Zustand)   │       │
│  └──────┬───────┘  └──────┬───────┘  └──────┬───────┘       │
│         │                  │                  │              │
│         └──────────────────┴──────────────────┘              │
│                            │                                 │
│                   ┌─────────▼─────────┐                      │
│                   │  engineAdapter    │                      │
│                   │  (lib/)           │                      │
│                   └─────────┬─────────┘                      │
└────────────────────────────┼─────────────────────────────────┘
                             │
                             │ @biosynth/engine
                             │
┌────────────────────────────▼─────────────────────────────────┐
│                      Engine Package                           │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐      │
│  │MoleculeGraph │  │LayoutEngine   │  │Molecule      │      │
│  │              │  │              │  │Serializer    │      │
│  └──────────────┘  └──────────────┘  └──────────────┘      │
└──────────────────────────────────────────────────────────────┘
                             │
                             │ HTTP API
                             │
┌────────────────────────────▼─────────────────────────────────┐
│                         Backend                               │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐       │
│  │   Routes     │  │   Models     │  │   Utils      │       │
│  │  /predict    │  │  Predictor   │  │  Featurizer  │       │
│  │  /generate   │  │              │  │              │       │
│  └──────────────┘  └──────────────┘  └──────────────┘       │
└──────────────────────────────────────────────────────────────┘
```

## Component Details

### Frontend (`frontend/`)

#### State Management
- **Store**: `src/store/moleculeStore.ts`
  - Manages current molecule, selected atoms, loading states
  - Handles API calls and predictions
  - Uses Zustand for global state

#### API Client
- **File**: `src/lib/api.ts`
  - Axios-based HTTP client
  - Functions: `predict()`, `generate()`, `predictFast()`
  - Reads API URL from `VITE_API_URL` env var

#### Engine Adapter
- **File**: `src/lib/engineAdapter.ts`
  - Converts `MoleculeGraph` to React-renderable format
  - Handles SMILES serialization via `MoleculeSerializer`
  - Bridges engine and frontend

#### Components
- **MoleculeViewer**: 3D visualization using React Three Fiber
- **Dashboard**: Main page with molecule viewer and controls

### Engine (`packages/engine/`)

#### Core Classes
- **MoleculeGraph**: Graph data structure for molecules
  - Atoms and bonds management
  - Formula and weight calculation
  - JSON serialization

- **LayoutEngine**: 3D geometry optimization
  - Force field calculations
  - Bond stretching forces
  - Non-bonded repulsion

- **MoleculeSerializer**: Serialization utilities
  - `toJSON()` / `fromJSON()`
  - `toSMILES()` / `fromSMILES()` (TODO: full implementation)

- **UndoStack**: History management
  - Undo/redo functionality
  - State snapshots

### Backend (`backend/`)

#### Structure
```
backend/
├── app.py              # FastAPI application
├── routes/             # API route handlers
│   ├── predict.py      # /predict endpoint
│   └── generate.py     # /generate endpoint
├── models/             # ML models
│   └── predictor.py    # PropertyPredictor
└── utils/              # Utilities
    └── featurizer.py   # SMILES featurization
```

#### API Endpoints

**POST /predict**
- Input: `{ smiles: string }`
- Output: `{ properties: { stability, toxicity, solubility, bioavailability, novelty } }`

**POST /generate**
- Input: `{ prompt: string }`
- Output: `{ smiles: string }`
- TODO: Implement transformer-based generation

**POST /predict-fast**
- Input: `{ smiles: string }`
- Output: `{ properties: {...} }`
- TODO: Implement ONNX inference

## Data Flow

### Prediction Flow
1. User modifies molecule in frontend
2. `MoleculeViewer` updates via Zustand store
3. `useEffect` in Dashboard triggers `fetchPredictions()`
4. Store calls `api.predict()` with SMILES
5. Backend `/predict` endpoint:
   - Validates SMILES
   - Featurizes using RDKit
   - Runs through PropertyPredictor
   - Returns properties
6. Store updates `backendPredictions`
7. Dashboard displays properties

### Generation Flow
1. User clicks "Generate Molecule"
2. Dashboard calls `generateMolecule()`
3. Store calls `api.generate()` with prompt
4. Backend `/generate` endpoint:
   - TODO: Use transformer model
   - Returns SMILES
5. Store creates MoleculeGraph from SMILES
6. MoleculeViewer renders new molecule

## Dependencies

### Frontend
- React 18.2.0
- TypeScript 5.3.3
- Zustand 4.5.0 (state)
- Axios 1.6.8 (HTTP)
- React Three Fiber 8.15.9 (3D)
- TailwindCSS 3.4.1 (styling)

### Engine
- TypeScript 5.3.3
- No external dependencies (pure TS)

### Backend
- FastAPI 0.110.0
- PyTorch 2.2.0
- RDKit 2023.09.2
- ONNX Runtime 1.17.0

## TODO Items

### High Priority
- [ ] Implement full SMILES serialization in `MoleculeSerializer`
- [ ] Implement SMILES parsing in `MoleculeSerializer`
- [ ] Replace DummyModel with real PyTorch model
- [ ] Add 3D interactions (atom selection, dragging)

### Medium Priority
- [ ] Implement transformer-based molecule generation
- [ ] Implement ONNX inference for fast predictions
- [ ] Add bond creation tool
- [ ] Add valence validation

### Low Priority
- [ ] Add comprehensive error handling
- [ ] Add loading states and spinners
- [ ] Add molecule library/explore page
- [ ] Add CI/CD pipeline

## 3D Interaction System

The molecular viewer implements a full 3D interaction system for manipulating molecules directly in the viewport.

### Architecture

```
┌─────────────────────────────────────────────────────────┐
│              MoleculeViewer Component                  │
│  ┌─────────────────────────────────────────────────┐   │
│  │         InteractionLayer                        │   │
│  │  - Drag handling                                │   │
│  │  - Raycasting                                   │   │
│  │  - Coordinate conversion                       │   │
│  └─────────────────────────────────────────────────┘   │
│  ┌─────────────────────────────────────────────────┐   │
│  │         SelectionManager (Singleton)            │   │
│  │  - hoveredAtomId                                │   │
│  │  - selectedAtomId                               │   │
│  │  - draggingAtomId                               │   │
│  │  - Event emitter pattern                        │   │
│  └─────────────────────────────────────────────────┘   │
│  ┌─────────────────────────────────────────────────┐   │
│  │         AtomMesh Components                     │   │
│  │  - Pointer event handlers                       │   │
│  │  - Visual feedback (hover/select/drag)         │   │
│  │  - Outlines for selection                      │   │
│  └─────────────────────────────────────────────────┘   │
│  ┌─────────────────────────────────────────────────┐   │
│  │         BondTool Hook                          │   │
│  │  - Two-atom selection logic                     │   │
│  │  - Bond creation                                │   │
│  │  - Escape key cancellation                     │   │
│  └─────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────┘
```

### Selection System

**SelectionManager** (`frontend/src/components/r3f/SelectionManager.ts`)
- Singleton pattern for global selection state
- Tracks three states: hover, select, drag
- Event emitter for reactive updates
- Communicates with Zustand store for persistence

**Selection Flow:**
1. User hovers over atom → `onHover(id)` → updates `hoveredAtomId`
2. User clicks atom → `onSelect(id)` → updates `selectedAtomId` + store
3. User drags atom → `startDrag(id)` → updates `draggingAtomId`

### Dragging System

**Raycasting** (`frontend/src/lib/raycasting.ts`)
- `getAtomUnderCursor()` - Finds atom mesh under mouse using raycasting
- `screenToWorld()` - Converts screen coordinates to 3D world coordinates
- Uses THREE.js Raycaster for intersection detection

**Drag Flow:**
1. User presses pointer down on atom → `startDrag(id)`
2. Pointer move events → `screenToWorld()` converts to 3D position
3. `updateAtomPosition()` updates MoleculeGraph
4. Store updates trigger re-render
5. Predictions automatically re-run on geometry change

**InteractionLayer Component:**
- Handles global drag state
- Manages invisible drag plane for coordinate conversion
- Listens to SelectionManager events
- Updates atom positions via `engineAdapter.updateAtomPosition()`

### Visual Feedback

**AtomMesh Visual States:**
- **Hover**: Blue outline glow (drei `<Outlines />`)
- **Selected**: Scale to 1.15x + outline
- **Dragging**: Opacity reduced to 0.7

**Cursor States:**
- Default: Normal cursor
- Hover: `pointer` cursor
- Dragging: `grabbing` cursor

### Bond Creation Tool

**BondTool Hook** (`frontend/src/components/r3f/BondTool.ts`)
- Monitors selected atom changes
- When two different atoms selected → creates bond
- Prevents duplicate bonds
- Escape key cancels bond creation

**Bond Creation Flow:**
1. User selects first atom → stored in `firstSelectedRef`
2. User selects second atom → `addBond(a1, a2)` called
3. MoleculeGraph updated → store updated → predictions re-run
4. UI updates to show new bond

### Rendering Pipeline

1. **MoleculeGraph** → `moleculeToRenderable()` → Renderable atoms/bonds
2. **AtomMesh** components render with interaction handlers
3. **SelectionManager** tracks state changes
4. **InteractionLayer** handles drag updates
5. **engineAdapter** syncs changes back to MoleculeGraph
6. **Store** triggers re-render and prediction updates

### Coordinate Systems

- **Screen Space**: Mouse coordinates (0-1 normalized)
- **World Space**: 3D coordinates in molecule space
- **Conversion**: Raycaster intersects with plane at atom's Y position

### Event Flow

```
User Action → AtomMesh Handler → SelectionManager → Store Update
                ↓
         InteractionLayer (drag)
                ↓
         engineAdapter.updateAtomPosition()
                ↓
         MoleculeGraph Update
                ↓
         Store.setMolecule()
                ↓
         Re-render + fetchPredictions()
```

### Performance Considerations

- Raycasting only on hover/click (not every frame)
- Position updates debounced via React state
- MoleculeGraph cloning only on actual changes
- Predictions re-run automatically but can be optimized with debouncing

### Testing

- **SelectionManager**: Unit tests for state tracking and events
- **Bond Management**: Engine tests for bond creation/removal
- **Position Updates**: Integration tests for drag-to-move

## ForceField Geometry Solver

The ForceField system provides molecular mechanics-based geometry optimization for realistic 3D conformations.

### Architecture

```
┌─────────────────────────────────────────────────────────┐
│              ForceField Engine                           │
│  ┌─────────────────────────────────────────────────┐     │
│  │         Energy Terms                            │     │
│  │  - Bond Stretching (E = k*(d-d0)²)            │     │
│  │  - Angle Bending (E = k*(θ-θ0)²)              │     │
│  │  - Lennard-Jones Repulsion                     │     │
│  └─────────────────────────────────────────────────┘     │
│  ┌─────────────────────────────────────────────────┐     │
│  │         Gradient Descent Optimizer              │     │
│  │  - Compute forces on all atoms                 │     │
│  │  - Update positions iteratively                │     │
│  │  - Minimize total energy                       │     │
│  └─────────────────────────────────────────────────┘     │
└─────────────────────────────────────────────────────────┘
```

### Energy Terms

**1. Bond Stretching**
```
E_stretch = k_bond * (d - d0)²
```
- `k_bond`: Force constant (kcal/mol/Å²)
- `d`: Current bond length
- `d0`: Equilibrium bond length (from covalent radii)
- Force: `F = -2k(d - d0)` along bond vector

**2. Angle Bending**
```
E_angle = k_angle * (θ - θ0)²
```
- `k_angle`: Force constant (kcal/mol/rad²)
- `θ`: Current angle
- `θ0`: Equilibrium angle (~109.5° for tetrahedral)
- Forces distributed perpendicular to bond vectors

**3. Lennard-Jones Repulsion**
```
E_lj = 4ε * ((σ/r)¹² - (σ/r)⁶)
```
- `ε`: Well depth (0.1 kcal/mol)
- `σ`: Van der Waals radius (3.0 Å)
- `r`: Distance between non-bonded atoms
- Only applied to non-bonded pairs

### Total Energy Equation

```
E_total = Σ E_stretch + Σ E_angle + Σ E_lj
```

### Gradient Descent Optimization

**Algorithm:**
1. Initialize forces for all atoms to zero
2. Compute bond stretching forces → add to atoms
3. Compute angle bending forces → add to atoms
4. Compute non-bonded repulsion forces → add to atoms
5. Update positions: `r_new = r_old + step * F`
6. Repeat for N iterations

**Parameters:**
- `iterations`: Number of optimization steps (default: 20)
- `step`: Learning rate / step size (default: 0.005)

### Integration Points

**After Atom Dragging:**
- `updateAtomPosition()` → quick relaxation (5 iterations)
- Smooths out distortions from manual dragging

**After Bond Creation:**
- `addBond()` → full relaxation (10 iterations)
- Ensures new bonds have proper geometry

**Future: WebWorker Optimization**
- TODO: Move optimization to WebWorker for non-blocking UI
- Currently runs asynchronously with `setTimeout`

### Vec3 Utility

**Vector Math Operations:**
- `add()`, `sub()`, `mul()`, `div()` - Basic arithmetic
- `dot()`, `cross()` - Vector products
- `length()`, `normalize()` - Magnitude operations
- `angleBetween()` - Angle calculation
- `distance()` - Distance between vectors

### Angle Detection

**MoleculeGraph Methods:**
- `getAnglesAround(atomId)`: Returns all angles centered on an atom
- `getAllAngles()`: Returns all unique angles in molecule
- Uses neighbor graph to find angle triples

### Performance Considerations

- Optimization runs asynchronously to avoid blocking UI
- Small iteration counts for real-time feedback (5-10 iterations)
- Full optimization can use more iterations (20-50) for final structure
- WebWorker implementation planned for large molecules

### Testing

**Unit Tests:**
- Bond force reduces stretched bonds
- Repulsion pushes atoms apart
- Optimization converges (energy decreases)
- Angle detection works correctly

## Development Guidelines

See `cursor.json` for detailed coding standards:
- Strong TypeScript typing (no `any`)
- Pure functions in engine
- Modular backend structure
- Aluminium design system
- Comprehensive testing

