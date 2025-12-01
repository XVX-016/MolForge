# Server Startup Guide

## Quick Start

1. **Stop any existing server**:
   ```powershell
   .\scripts\stop_server.ps1
   # Or manually:
   Get-NetTCPConnection -LocalPort 8000 -ErrorAction SilentlyContinue | Select-Object -ExpandProperty OwningProcess -Unique | ForEach-Object { Stop-Process -Id $_ -Force }
   ```

2. **Start the server**:
   ```powershell
   cd backend
   $env:PYTHONPATH = "C:\Computing\MolForge"
   python -m uvicorn app:app --host 0.0.0.0 --port 8000
   ```

3. **Test the API** (in another terminal):
   ```powershell
   cd backend
   python scripts/test_api_endpoints.py
   ```

## What Was Fixed

### 1. PyTorch DLL Loading Error
- **Problem**: PyTorch DLL fails to load at module import time
- **Solution**: Implemented lazy imports - PyTorch is only loaded when actually needed
- **Files**: `ml/featurize.py`, `ml/prediction_engine.py`, `ml/registry.py`

### 2. Pydantic Warnings
- **Problem**: `Field "model_id" has conflict with protected namespace "model_"`
- **Solution**: Added `Config` class with `protected_namespaces = ()`
- **Files**: `api/predict.py`

### 3. Port Already in Use
- **Problem**: Port 8000 already bound
- **Solution**: Created `scripts/stop_server.ps1` to kill processes on port 8000

## Expected Behavior

When server starts, you should see:
- ✅ No PyTorch DLL errors (lazy import)
- ✅ No Pydantic warnings
- ✅ Server starts on port 8000
- ✅ API endpoints respond correctly

## If Errors Persist

See `TROUBLESHOOTING.md` for detailed solutions.

