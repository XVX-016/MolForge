# Troubleshooting Guide

## PyTorch DLL Loading Error

**Error**: `[WinError 1114] A dynamic link library (DLL) initialization routine failed`

**Solution**: 
- PyTorch is installed but DLL loading fails at import time
- Fixed with lazy imports - PyTorch is now loaded only when needed
- If error persists, try:
  1. Reinstall PyTorch: `pip uninstall torch torchvision torchaudio && pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu`
  2. Install Visual C++ Redistributables
  3. Restart Python/terminal

## Port Already in Use

**Error**: `[Errno 10048] error while attempting to bind on address ('0.0.0.0', 8000)`

**Solution**:
```powershell
# Stop process on port 8000
Get-NetTCPConnection -LocalPort 8000 -ErrorAction SilentlyContinue | Select-Object -ExpandProperty OwningProcess -Unique | ForEach-Object { Stop-Process -Id $_ -Force }

# Or use the script
.\backend\scripts\stop_server.ps1
```

## Test Script Path Issues

**Error**: `can't open file 'C:\\Computing\\MolForge\\backend\\backend\\scripts\\...'`

**Solution**: 
- If in `backend/` directory, use: `python scripts/test_api_endpoints.py`
- If in root, use: `python backend/scripts/test_api_endpoints.py`

## Pydantic Warnings

**Warning**: `Field "model_id" has conflict with protected namespace "model_"`

**Solution**: Fixed by adding `Config` class with `protected_namespaces = ()`

## Server Startup

**Correct startup command**:
```powershell
cd backend
$env:PYTHONPATH = "C:\Computing\MolForge"
python -m uvicorn app:app --host 0.0.0.0 --port 8000
```

