# PowerShell script to start the backend server
# Usage: .\start_backend.ps1
# Run from the backend/ directory

Write-Host "Starting MolForge Backend Server..." -ForegroundColor Green
Write-Host ""

# Check if we're in the backend directory
if (-not (Test-Path "app.py")) {
    Write-Host "Error: app.py not found. Please run this script from the backend/ directory" -ForegroundColor Red
    exit 1
}

# Activate virtual environment (check both .venv and ../.venv)
if (Test-Path ".venv\Scripts\Activate.ps1") {
    .\.venv\Scripts\Activate.ps1
    Write-Host "Virtual environment activated (.venv)" -ForegroundColor Green
} elseif (Test-Path "..\.venv\Scripts\Activate.ps1") {
    ..\.venv\Scripts\Activate.ps1
    Write-Host "Virtual environment activated (../.venv)" -ForegroundColor Green} else {
    Write-Host "No virtual environment found. Using system Python..." -ForegroundColor Yellow
    Write-Host "   Tip: Create one with: python -m venv .venv" -ForegroundColor Yellow
}

# Check if uvicorn is installed
python -c "import uvicorn" 2>$null
if ($LASTEXITCODE -ne 0) {
    Write-Host "uvicorn not installed. Installing..." -ForegroundColor Yellow
    pip install uvicorn
}

# Check if .env exists, create from example if not
if (-not (Test-Path ".env")) {
    if (Test-Path "env.example") {
        Write-Host "Creating .env from env.example..." -ForegroundColor Yellow
        Copy-Item "env.example" ".env"
    }
}

# Set PYTHONPATH to include parent directory so backend imports work
# Get the parent directory (repo root)
$currentDir = Get-Location
$parentDir = Split-Path -Parent $currentDir
$env:PYTHONPATH = "$parentDir;$env:PYTHONPATH"

Write-Host "PYTHONPATH set to: $parentDir" -ForegroundColor Gray

# Start the server
Write-Host ""
Write-Host "Starting server on http://localhost:8000" -ForegroundColor Cyan
Write-Host "API docs available at http://localhost:8000/docs" -ForegroundColor Cyan
Write-Host "Health check: http://localhost:8000/health" -ForegroundColor Cyan
Write-Host ""
Write-Host "Press Ctrl+C to stop the server" -ForegroundColor Yellow
Write-Host ""

# Run from backend directory but with PYTHONPATH set to parent
# Use --reload-dir to only watch current directory and avoid path issues
python -m uvicorn app:app --reload --reload-dir . --host 127.0.0.1 --port 8000

