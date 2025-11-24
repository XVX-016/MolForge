# ML Architecture & Logic Flow Design

## 1️⃣ Overall Logic Flow

### Input Layer
- User uploads/draws molecule → Parsed into graph representation (atoms + bonds)
- Molecular descriptors computed immediately

### Preprocessing & Feature Extraction
- **Structural features:** atom types, bond types, rings, functional groups
- **Spectroscopy features:** IR active groups, NMR shielding predictions
- **Energy features:** MM energies, bond lengths, angles
- **Quantum features:** HOMO/LUMO, ESP maps

### Engine Calculations (Parallel)
- Spectroscopy Engines: IR, NMR, MassSpec
- KAB Analysis: Binding sites, toxicity alerts, stability rules
- Energy & Geometry: Minimized geometry, molecular dynamics
- Quantum Calculations: HOMO/LUMO, ESP

### ML Prediction Layer
- Input: Feature vector from preprocessing
- Output: Predicted properties (solubility, toxicity, bioactivity, etc.)
- Per-atom/per-bond contributions for explainability

### Integration & Visualization
- Unified API aggregates all engine outputs
- Frontend visualizes: 3D molecule, highlights, spectra, predictions

### Optional Feedback Loop
- ML predictions can flag KAB alerts or adjust geometry

---

## 2️⃣ Feature Engineering for ML

### Graph-based Features
- **Node features:** atom type, hybridization, charge, aromaticity, electronegativity
- **Edge features:** bond type, conjugation, ring membership

### Chemical Descriptors
- Molecular weight, logP, TPSA, H-bond donors/acceptors
- Ring counts, rotatable bonds, functional group counts

### Engine-derived Features
- IR peak counts per functional group
- NMR shielding averages
- MM energy components (bond, angle, torsion)
- HOMO/LUMO gap

### Optional Embeddings
- Pre-trained molecular embeddings (ChemBERTa, GNN embeddings)

---

## 3️⃣ ML Model Design

### Model Type
- **Graph Neural Network (GNN)** for per-atom/bond predictions
- **Feedforward Neural Network (FNN)** for global property predictions

### Input/Output
- **Input:** graph + features vector
- **Output:** property predictions + per-atom contributions

### Explainability
- **Integrated Gradients** or **attention weights** for per-atom scores
- Used to highlight atoms/bonds in 3D viewer

### Training
- Dataset: PubChem, QM9, or custom chemistry datasets
- Targets: solubility, toxicity, reactivity, synthetic accessibility
- Loss: MSE for regression, cross-entropy for classification

---

## 4️⃣ Execution Logic

```
User uploads molecule → Parser → Feature extractor → 
[Engines: IR, NMR, MassSpec, KAB, Energy, Geometry, Quantum] → 
ML feature vector → ML prediction → Explanation engine → 
Unified API → Frontend 3D visualization + spectra + predictions + highlights
```

### Notes
- Engines run in parallel where possible
- ML depends on feature extraction but can run asynchronously
- Use caching for repeated calculations
- ML predictions trigger visual cues (red=unstable, yellow=moderate toxicity)

---

## 5️⃣ Implementation Structure

### Backend
- `backend/chem/ml/features.py` - Feature extraction
- `backend/chem/ml/predict.py` - ML model inference
- `backend/chem/ml/explain.py` - Explainability engine
- `backend/chem/ml/models/` - Trained model files (optional)

### Frontend
- `frontend/src/components/ml/MLGraph.tsx` - Visualization
- `frontend/src/components/ml/PredictionExplanation.tsx` - Explanations
- `frontend/src/components/lab/panels/MLPanel.tsx` - Lab integration

### API
- `POST /api/ml/predict` - Property predictions
- `POST /api/ml/features` - Feature extraction
- `POST /api/ml/explain` - AI explanations

