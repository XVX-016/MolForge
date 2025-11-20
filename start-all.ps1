# Simple script to start both backend and frontend
# Run from repo root: .\start-all.ps1

Write-Host "🚀 Starting Backend and Frontend..." -ForegroundColor Green
Write-Host ""

# Start backend in background
Write-Host "🔧 Starting Backend (new window)..." -ForegroundColor Cyan
Start-Process powershell -ArgumentList "-NoExit", "-Command", "cd '$PSScriptRoot'; .\start-backend.ps1"

# Wait a bit
Start-Sleep -Seconds 2

# Start frontend in current window
Write-Host "🎨 Starting Frontend..." -ForegroundColor Cyan
Write-Host ""

# Install dependencies if needed
if (-not (Test-Path "frontend\node_modules")) {
    Write-Host "📦 Installing frontend dependencies..." -ForegroundColor Yellow
    Set-Location frontend
    npm install
    Set-Location ..
}

Write-Host ""
Write-Host "✅ Servers starting!" -ForegroundColor Green
Write-Host "📍 Backend:  http://localhost:8000" -ForegroundColor Cyan
Write-Host "📍 Frontend: http://localhost:5173" -ForegroundColor Cyan
Write-Host ""
Write-Host "Press Ctrl+C to stop frontend" -ForegroundColor Yellow
Write-Host ""

Set-Location frontend
npm run dev

