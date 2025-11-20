# PowerShell script to start the backend server
# Usage: .\start-backend.ps1
# Run from the repo root

Write-Host "🚀 Starting MolForge Backend Server..." -ForegroundColor Green
Write-Host ""

# Check if we're in the repo root
if (-not (Test-Path "backend\app.py")) {
    Write-Host "❌ Error: backend\app.py not found. Please run this script from the repo root" -ForegroundColor Red
    exit 1
}

# Activate virtual environment (check both backend/.venv and .venv)
if (Test-Path "backend\.venv\Scripts\Activate.ps1") {
    . backend\.venv\Scripts\Activate.ps1
    Write-Host "✅ Virtual environment activated (backend/.venv)" -ForegroundColor Green
} elseif (Test-Path ".venv\Scripts\Activate.ps1") {
    . .venv\Scripts\Activate.ps1
    Write-Host "✅ Virtual environment activated (.venv)" -ForegroundColor Green
} else {
    Write-Host "⚠️  No virtual environment found. Using system Python..." -ForegroundColor Yellow
    Write-Host "   Tip: Create one with: cd backend; python -m venv .venv" -ForegroundColor Yellow
}

# Check if uvicorn is installed
python -c "import uvicorn" 2>$null
if ($LASTEXITCODE -ne 0) {
    Write-Host "❌ uvicorn not installed. Installing..." -ForegroundColor Yellow
    pip install uvicorn
}

# Check if .env exists, create from example if not
if (-not (Test-Path "backend\.env")) {
    if (Test-Path "backend\env.example") {
        Write-Host "📝 Creating backend/.env from env.example..." -ForegroundColor Yellow
        Copy-Item "backend\env.example" "backend\.env"
    }
}

# Start the server from repo root (so backend imports work)
Write-Host ""
Write-Host "🚀 Starting server on http://localhost:8000" -ForegroundColor Cyan
Write-Host "📚 API docs available at http://localhost:8000/docs" -ForegroundColor Cyan
Write-Host "❤️  Health check: http://localhost:8000/health" -ForegroundColor Cyan
Write-Host ""
Write-Host "Press Ctrl+C to stop the server" -ForegroundColor Yellow
Write-Host ""

# Run from repo root with backend.app:app
python -m uvicorn backend.app:app --reload --host 0.0.0.0 --port 8000

