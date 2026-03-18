Write-Host "[INFO] Starting MolForge Backend Server..." -ForegroundColor Green
Write-Host ""

if (-not (Test-Path "backend\app.py")) {
    Write-Host "[ERROR] backend\app.py not found. Please run this script from the repo root." -ForegroundColor Red
    exit 1
}

function Test-PythonInterpreter {
    param([string]$Command)

    try {
        & $Command --version 1>$null 2>$null
        return $LASTEXITCODE -eq 0
    } catch {
        return $false
    }
}

$repoRoot = (Get-Location).Path
$backendPath = Join-Path $repoRoot "backend"
$candidateInterpreters = @(
    @{ Path = (Join-Path $backendPath ".venv\Scripts\python.exe"); Label = "backend/.venv" },
    @{ Path = (Join-Path $repoRoot ".venv\Scripts\python.exe"); Label = ".venv" },
    @{ Path = "python"; Label = "system Python" }
)

$pythonExe = $null
foreach ($candidate in $candidateInterpreters) {
    $candidatePath = $candidate.Path
    if ($candidatePath -ne "python" -and -not (Test-Path $candidatePath)) {
        continue
    }

    if (Test-PythonInterpreter -Command $candidatePath) {
        $pythonExe = $candidatePath
        Write-Host "[OK] Using $($candidate.Label): $candidatePath" -ForegroundColor Green
        break
    }

    if ($candidatePath -ne "python") {
        Write-Host "[WARN] Skipping broken interpreter at $($candidate.Label): $candidatePath" -ForegroundColor Yellow
    }
}

if (-not $pythonExe) {
    Write-Host "[ERROR] No usable Python interpreter was found." -ForegroundColor Red
    Write-Host "        Recreate the backend venv with Python 3.10, for example:" -ForegroundColor Red
    Write-Host "        py -3.10 -m venv backend\.venv" -ForegroundColor Red
    exit 1
}

$pythonVersion = & $pythonExe -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}')" 2>$null
if (-not $pythonVersion) {
    Write-Host "[ERROR] Failed to query the selected Python interpreter: $pythonExe" -ForegroundColor Red
    exit 1
}

if (-not $pythonVersion.StartsWith("3.10")) {
    Write-Host "[WARN] Active interpreter is Python $pythonVersion. This repo is now aligned to Python 3.10." -ForegroundColor Yellow
}

if ((Test-Path "backend\.venv\Scripts\Activate.ps1") -and -not (Test-Path "backend\.venv\Scripts\python.exe")) {
    Write-Host "[WARN] backend/.venv still looks like an old 3.11 environment shell. Recreate it if imports stay flaky." -ForegroundColor Yellow
}

if (-not (Test-Path "backend\.env")) {
    if (Test-Path "backend\env.example") {
        Write-Host "[INFO] Creating backend/.env from env.example..." -ForegroundColor Yellow
        Copy-Item "backend\env.example" "backend\.env"
    } else {
        Write-Host "[WARN] backend/env.example missing; skipping .env creation." -ForegroundColor Yellow
    }
}

$env:PYTHONPATH = "$repoRoot;$backendPath"

$missingModules = & $pythonExe -c "import importlib.util as u; mods=['uvicorn','fastapi','sqlmodel','dotenv']; missing=[m for m in mods if u.find_spec(m) is None]; print(','.join(missing))"
if ($missingModules) {
    Write-Host "[WARN] Missing backend dependencies: $missingModules" -ForegroundColor Yellow
    if (Test-Path "backend\requirements.txt") {
        Write-Host "[INFO] Installing backend requirements into the selected interpreter..." -ForegroundColor Yellow
        & $pythonExe -m pip install -r "backend\requirements.txt"
        if ($LASTEXITCODE -ne 0) {
            Write-Host "[ERROR] Dependency installation failed. Please fix the environment and rerun this script." -ForegroundColor Red
            exit $LASTEXITCODE
        }
    } else {
        Write-Host "[ERROR] backend/requirements.txt not found, cannot install dependencies." -ForegroundColor Red
        exit 1
    }
}

Write-Host ""
Write-Host "[INFO] Starting server on http://localhost:8000" -ForegroundColor Cyan
Write-Host "[INFO] API docs available at http://localhost:8000/docs" -ForegroundColor Cyan
Write-Host "[INFO] Health check: http://localhost:8000/health" -ForegroundColor Cyan
Write-Host ""
Write-Host "Press Ctrl+C to stop the server" -ForegroundColor Yellow
Write-Host ""

& $pythonExe -m uvicorn backend.app:app --reload --reload-dir backend --host 127.0.0.1 --port 8000
