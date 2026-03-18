# MolForge

MolForge is a monorepo for molecular design and chemistry workflows. It combines a React/Vite frontend, a FastAPI backend, and a shared TypeScript engine package used by the UI.

## What Is In This Repo

- `frontend/`: React 18 + Vite application for dashboard, library, studio, and lab flows.
- `backend/`: FastAPI service for molecule, search, spectroscopy, energy, reaction, retrosynthesis, screening, dashboard, and studio APIs.
- `packages/engine/`: shared TypeScript package consumed as `@biosynth/engine`.
- `data/`, `models/`, `uploads/`: local project assets and generated artifacts.
- `biosynth.db`: checked-in SQLite database used for local development.

## Tech Stack

- Frontend: React, TypeScript, Vite, Zustand, Three.js, Vitest, ESLint
- Backend: Python 3.10, FastAPI, SQLModel, Uvicorn, pytest
- Chemistry/ML: RDKit, NumPy, ONNX Runtime, optional model/tooling dependencies

## Strict Local Requirements

- Node.js 20+ recommended
- npm 10+ recommended
- Python 3.10 required for the backend runtime
- PowerShell for the provided startup script on Windows
- RDKit is required for chemistry-heavy backend endpoints

Do not use Python 3.11 for the local backend environment unless you intentionally re-align the project. The repository has been normalized for Python 3.10.

## Quick Start

### 1. Install frontend and workspace packages

```powershell
npm install
```

### 2. Create a backend virtual environment

```powershell
py -3.10 -m venv backend\.venv
backend\.venv\Scripts\python.exe -m pip install --upgrade pip
backend\.venv\Scripts\python.exe -m pip install -r backend\requirements.txt
```

If you need the larger ML stack as well:

```powershell
backend\.venv\Scripts\python.exe -m pip install -r backend\requirements-model.txt
```

### 3. Start the backend

Preferred:

```powershell
powershell -ExecutionPolicy Bypass -File .\start-backend.ps1
```

Alternative workspace command:

```powershell
npm --workspace backend run dev
```

### 4. Start the frontend

```powershell
npm --workspace frontend run dev
```

### 5. Start the full stack

```powershell
npm run dev
```

## Key URLs

- Frontend dev server: `http://127.0.0.1:5173`
- Backend API: `http://127.0.0.1:8000`
- OpenAPI docs: `http://127.0.0.1:8000/docs`
- Health check: `http://127.0.0.1:8000/health`

## Project Commands

```powershell
npm run dev
npm run build
npm test
npm --workspace frontend run dev
npm --workspace frontend run build
npm --workspace frontend run test
npm --workspace frontend run lint
npm --workspace backend run dev
python -m pytest -q
python -m pytest -q backend/tests
```

## Backend Runtime Notes

- The backend startup script now prefers a valid interpreter and skips stale or broken virtual environments.
- `backend/app.py` loads routers defensively so missing optional third-party dependencies do not always prevent the entire API from booting.
- Some chemistry routes still need RDKit and related packages at runtime. Missing optional dependencies will disable those routes rather than silently emulating chemistry behavior.

## Testing Guidance

- Frontend TypeScript currently compiles with `tsc`.
- Frontend ESLint has a large existing backlog across many files. Treat lint cleanup as incremental work, not a single-step fix.
- Pytest is configured to focus on `backend/tests` and ignore `backend/.venv` and other generated folders.

## Git Hygiene

- Do not commit local virtual environments, caches, weights, or generated datasets.
- Keep large local-only assets out of Git unless they are intentionally versioned.
- This repository may contain user-owned uncommitted deployment file deletions. Do not blindly restore them unless that is the intended change.

## Known Hotspots

- Backend duplication risk: `backend/routes/` and `backend/api/`
- Frontend duplication risk: `frontend/src/lib/molecule`, `frontend/src/chemcore`, `frontend/src/kernel`, `frontend/src/utils`, and engine-related code
- `frontend/vite.config.ts` resolves `@biosynth/engine` from built output, so engine changes may require rebuilding the package

## Recommended Cleanup Order

1. Stabilize backend environment on Python 3.10.
2. Keep startup scripts and workspace commands aligned.
3. Tackle frontend lint debt by area, not file-by-file at random.
4. Consolidate overlapping molecule models before large feature work.

