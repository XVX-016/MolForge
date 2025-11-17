# Backend Structure

This document describes the structure and organization of the BioSynth AI backend.

## Directory Structure

```
backend/
├── app.py                 # FastAPI application entry point
├── db.py                  # Database configuration and session management
├── models_db.py           # SQLModel database models
├── requirements.txt       # Python dependencies
├── Dockerfile            # Docker container configuration
├── run_backend.sh        # Shell script to run the backend
│
├── models/               # ML models and predictors
│   ├── __init__.py
│   ├── predictor.py      # Main property predictor
│   ├── property_predictor.py  # PyTorch neural network model
│   ├── onnx_predictor.py # ONNX-based fast inference
│   ├── export_to_onnx.py # Export PyTorch model to ONNX
│   └── train_property_predictor.py  # Training script
│
├── routes/               # API route handlers
│   ├── __init__.py
│   ├── predict.py        # Property prediction endpoints
│   ├── generate.py       # Molecule generation endpoints
│   ├── library.py        # Molecule library management
│   └── admin.py          # Admin endpoints
│
├── utils/                # Utility functions
│   ├── __init__.py
│   └── featurizer.py     # SMILES featurization (RDKit)
│
├── jobs/                 # Background job processing
│   └── worker.py         # RQ worker for async tasks
│
├── simulation/           # Reaction simulation
│   └── reaction_engine.py  # Chemical reaction engine
│
└── test_*.py            # Test files
```

## Key Components

### 1. Application Entry Point (`app.py`)
- FastAPI application setup
- CORS middleware configuration
- Router mounting
- Health check endpoint

### 2. Database (`db.py`, `models_db.py`)
- SQLModel ORM setup
- Database connection management
- Molecule model definition

### 3. ML Models (`models/`)
- **PropertyPredictor**: PyTorch neural network for property prediction
- **ONNXPredictor**: Optimized ONNX inference
- Training and export utilities

### 4. API Routes (`routes/`)
- `/predict`: Property prediction endpoints
- `/generate`: Molecule generation
- `/molecules`: Library management
- `/admin`: Administrative functions

### 5. Utilities (`utils/`)
- SMILES string validation and featurization
- Morgan fingerprint generation using RDKit

### 6. Background Jobs (`jobs/`)
- Redis Queue (RQ) worker setup
- Async task processing for:
  - Molecule generation
  - Property prediction
  - Property optimization
  - Reaction processing

## Environment Variables

See `.env.example` for required environment variables:
- `DATABASE_URL`: Database connection string
- `REDIS_URL`: Redis connection for job queue
- `ENVIRONMENT`: Application environment (development/production)

## Dependencies

### Core Framework
- FastAPI: Web framework
- Uvicorn: ASGI server
- Pydantic: Data validation
- SQLModel: Database ORM

### ML/AI
- PyTorch: Deep learning framework
- Transformers: Hugging Face transformers
- NumPy: Numerical computing
- ONNXRuntime: Optimized inference

### Chemistry
- RDKit: Chemical informatics toolkit
- Pandas: Data manipulation

### Background Jobs
- Redis: Job queue backend
- RQ: Python job queue library

## Running the Backend

### Development
```bash
# Create virtual environment
python3.11 -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Set up environment variables
cp .env.example .env
# Edit .env with your configuration

# Run the server
uvicorn app:app --reload --host 0.0.0.0 --port 8000
```

### Docker
```bash
docker build -t biosynth-backend .
docker run -p 8000:8000 --env-file .env biosynth-backend
```

### Background Worker
```bash
# Start Redis (if not running)
redis-server

# Start worker
python -m backend.jobs.worker
```

## API Endpoints

- `GET /`: API information
- `GET /health`: Health check
- `POST /predict-fast`: Fast property prediction (ONNX)
- `/predict/*`: Property prediction routes
- `/generate/*`: Molecule generation routes
- `/molecules/*`: Library management routes
- `/admin/*`: Admin routes

## Testing

Run tests with:
```bash
pytest
```

Individual test files:
- `test_featurizer.py`: SMILES featurization tests
- `test_predict_*.py`: Prediction endpoint tests
- `test_library.py`: Library management tests

