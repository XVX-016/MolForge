# Backend ML System Documentation

## Overview

The backend ML system provides molecular property prediction using deep neural networks. It includes training pipelines, model inference, and ONNX-accelerated prediction endpoints.

## Model Architecture

### PropertyPredictor Neural Network

**Architecture:**
```
Input: 2048-bit Morgan fingerprint
  ↓
FC1: 2048 → 1024 (ReLU + Dropout 0.3)
  ↓
FC2: 1024 → 512 (ReLU + Dropout 0.3)
  ↓
FC3: 512 → 128 (ReLU + Dropout 0.3)
  ↓
Output: 128 → 5 properties
```

**Output Properties:**
1. `stability` - Molecular stability score (0-1)
2. `toxicity` - Toxicity score (0-1)
3. `solubility` - Solubility score (0-1)
4. `bioavailability` - Bioavailability score (0-1)
5. `novelty` - Novelty score (0-1)

**Key Features:**
- Kaiming weight initialization
- Dropout regularization (0.3)
- ReLU activations
- No output activation (regression)

## Dataset Format

### Training CSV Structure

**File:** `backend/data/molecules.csv`

**Columns:**
- `smiles` (string): SMILES representation of molecule
- `stability` (float): Stability score (0-1)
- `toxicity` (float): Toxicity score (0-1)
- `solubility` (float): Solubility score (0-1)
- `bioavailability` (float): Bioavailability score (0-1)
- `novelty` (float): Novelty score (0-1)

**Example:**
```csv
smiles,stability,toxicity,solubility,bioavailability,novelty
CCO,0.75,0.25,0.80,0.70,0.10
CC,0.85,0.15,0.60,0.65,0.05
c1ccccc1,0.70,0.30,0.50,0.60,0.20
```

**Dummy Dataset Generation:**
If the dataset file doesn't exist, `train_property_predictor.py` will automatically generate a dummy dataset with 50 samples for testing purposes.

## Training Pipeline

### Training Script

**File:** `backend/models/train_property_predictor.py`

**Usage:**
```bash
python -m backend.models.train_property_predictor \
    --csv backend/data/molecules.csv \
    --epochs 50 \
    --batch-size 8 \
    --lr 0.001 \
    --output backend/weights/property_predictor.pt
```

**Training Process:**
1. Load dataset from CSV
2. Generate Morgan fingerprints (2048-bit) for each SMILES
3. Create DataLoader with batching
4. Train for N epochs with:
   - Adam optimizer
   - MSELoss (regression)
   - Learning rate scheduling
5. Save model weights to `backend/weights/property_predictor.pt`

**Training Parameters:**
- **Epochs:** 50 (default)
- **Batch Size:** 8 (default)
- **Learning Rate:** 0.001 (default)
- **Loss Function:** MSELoss
- **Optimizer:** Adam

## Fingerprint Generation

### Morgan Fingerprints

**Implementation:** `backend/utils/featurizer.py`

**Method:** `featurize_smiles(smiles: str) -> List[float]`

**Parameters:**
- **Radius:** 2 (bond distance)
- **Bits:** 2048
- **Algorithm:** RDKit Morgan fingerprint

**Process:**
1. Parse SMILES string using RDKit
2. Generate Morgan fingerprint with radius=2
3. Convert to bit vector (2048 bits)
4. Return as list of floats (0.0 or 1.0)

## Model Inference

### PyTorch Inference

**File:** `backend/models/predictor.py`

**Class:** `ModelLoader`

**Usage:**
```python
from backend.models.predictor import get_predictor

predictor = get_predictor()
features = featurize_smiles("CCO")
properties = predictor.predict(features)
```

**Process:**
1. Load model weights (or use untrained model if weights missing)
2. Convert features to PyTorch tensor
3. Run forward pass (no gradient computation)
4. Convert output to Python dict

### ONNX Inference

**File:** `backend/models/onnx_predictor.py`

**Class:** `ONNXPredictor`

**Export:**
```bash
python -m backend.models.export_to_onnx \
    --weights backend/weights/property_predictor.pt \
    --output backend/weights/predictor.onnx
```

**Usage:**
```python
from backend.models.onnx_predictor import get_onnx_predictor

predictor = get_onnx_predictor()
features = featurize_smiles("CCO")
properties = predictor.predict(features)
```

**Benefits:**
- Faster inference (optimized runtime)
- Smaller model size
- Cross-platform compatibility
- No PyTorch dependency at inference time

## API Endpoints

### POST /predict

**Request:**
```json
{
  "smiles": "CCO"
}
```

**Response:**
```json
{
  "properties": {
    "stability": 0.75,
    "toxicity": 0.25,
    "solubility": 0.80,
    "bioavailability": 0.70,
    "novelty": 0.10
  }
}
```

### POST /predict-fast

**Request:**
```json
{
  "smiles": "CCO"
}
```

**Response:**
```json
{
  "properties": {
    "stability": 0.75,
    "toxicity": 0.25,
    "solubility": 0.80,
    "bioavailability": 0.70,
    "novelty": 0.10
  }
}
```

**Note:** Falls back to PyTorch predictor if ONNX model not available.

### GET /predict/last

**Response:**
```json
{
  "properties": {
    "stability": 0.75,
    "toxicity": 0.25,
    "solubility": 0.80,
    "bioavailability": 0.70,
    "novelty": 0.10
  }
}
```

**Use Case:** Debugging - returns last predictions made.

### POST /predict/bulk

**Request:**
```json
{
  "smiles_list": ["CCO", "CC", "C"]
}
```

**Response:**
```json
{
  "results": [
    {
      "smiles": "CCO",
      "error": null,
      "properties": {
        "stability": 0.75,
        "toxicity": 0.25,
        "solubility": 0.80,
        "bioavailability": 0.70,
        "novelty": 0.10
      }
    },
    {
      "smiles": "CC",
      "error": null,
      "properties": {...}
    },
    {
      "smiles": "C",
      "error": null,
      "properties": {...}
    }
  ]
}
```

**Use Case:** Batch prediction for multiple molecules.

## ONNX Acceleration

### Export Process

**File:** `backend/models/export_to_onnx.py`

**Features:**
- Dynamic batch size support
- Opset version 11
- Constant folding optimization
- Model verification
- Test inference after export

**Dynamic Axes:**
- `fingerprint`: batch dimension (0)
- `properties`: batch dimension (0)

### ONNX Runtime

**Benefits:**
- 2-5x faster inference than PyTorch
- Lower memory footprint
- Optimized for CPU inference
- Cross-platform (Windows, Linux, macOS)

**Performance:**
- PyTorch: ~5-10ms per prediction
- ONNX: ~2-5ms per prediction

## Testing

### Unit Tests

**Files:**
- `backend/test_predictor_load.py` - Model loading tests
- `backend/test_predict_endpoint.py` - API endpoint tests
- `backend/test_onnx_predict.py` - ONNX inference tests

**Test Coverage:**
- Model initialization
- Forward pass correctness
- Output shape validation
- API endpoint responses
- ONNX vs PyTorch consistency

### Running Tests

```bash
# All ML tests
pytest backend/test_predictor_load.py backend/test_predict_endpoint.py backend/test_onnx_predict.py

# Individual test files
pytest backend/test_predictor_load.py -v
```

## Workflow

### Complete Training & Deployment Workflow

1. **Prepare Dataset:**
   ```bash
   # Create CSV with molecules and properties
   # Or use auto-generated dummy dataset
   ```

2. **Train Model:**
   ```bash
   python -m backend.models.train_property_predictor \
       --csv backend/data/molecules.csv \
       --epochs 50
   ```

3. **Export to ONNX:**
   ```bash
   python -m backend.models.export_to_onnx \
       --weights backend/weights/property_predictor.pt \
       --output backend/weights/predictor.onnx
   ```

4. **Start Backend:**
   ```bash
   cd backend
   uvicorn app:app --reload
   ```

5. **Test Endpoints:**
   ```bash
   curl -X POST http://localhost:8000/predict \
       -H "Content-Type: application/json" \
       -d '{"smiles": "CCO"}'
   ```

## File Structure

```
backend/
├── models/
│   ├── property_predictor.py      # PyTorch model definition
│   ├── train_property_predictor.py # Training script
│   ├── predictor.py               # ModelLoader & PropertyPredictor
│   ├── export_to_onnx.py          # ONNX export script
│   └── onnx_predictor.py           # ONNX inference
├── routes/
│   └── predict.py                  # API endpoints
├── utils/
│   └── featurizer.py               # SMILES → fingerprint
├── weights/
│   ├── property_predictor.pt       # PyTorch weights
│   └── predictor.onnx              # ONNX model
└── data/
    └── molecules.csv                # Training dataset
```

## Future Enhancements

- [ ] Multi-task learning with shared encoder
- [ ] Graph neural networks (GNN) for better molecular representation
- [ ] Transfer learning from pre-trained models
- [ ] Hyperparameter tuning with Optuna
- [ ] Model versioning and A/B testing
- [ ] Real-time model updates without restart
- [ ] GPU acceleration for training
- [ ] Distributed training for large datasets

