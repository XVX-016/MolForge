# MolForge Studio Contract

MolForge Studio is an AI-orchestrated molecular design environment.
All molecule mutations must comply with this contract.

---

## Molecule Graph Schema

### Atom
- **id**: string (UUID)
- **element**: string (e.g., "C", "N", "O", "H")
- **valence**: number (standard chemical valence)
- **position**: `[number, number, number]` (3D coordinates)

### Bond
- **id**: string (UUID)
- **from**: Atom.id
- **to**: Atom.id
- **order**: `1 | 2 | 3 | "aromatic"`

### Molecule
- **atoms**: Atom[]
- **bonds**: Bond[]

---

## Studio Modes

### DESIGN
- **Allowed**: Topology changes (add/remove/replace atoms and bonds).
- **Disallowed**: Real-time geometry optimization (handled separately).

### OPTIMIZE
- **Allowed**: Atom position updates only.
- **Disallowed**: Adding or removing atoms/bonds.

### SIMULATE
- **Read-only**: Structure is locked.
- **Output**: Generates trajectory data or reaction pathways.

---

## Action Flow

1. **User Input**: Natural language command in `AIActionCard`.
2. **Gemini Processing**: Intent parsing -> Action JSON generation.
3. **Validation Gate**: `MutationValidator` checks action against the Contract and Mode.
4. **State Mutation**: `studioStore` updates the `moleculeGraph`.
5. **System Feedback**: Results and warnings are logged in `AILogCard`.

---

## Failure Handling

- **Invalid Action**: Rejection with a clear chemical explanation in the log.
- **Hallucination**: Actions not in the `StudioAction` vocabulary are ignored.
- **Constraint Violation**: Unauthorized mutations (e.g., adding atoms in OPTIMIZE mode) are blocked.
