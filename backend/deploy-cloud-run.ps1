# PowerShell script to deploy to Google Cloud Run
# Usage: .\deploy-cloud-run.ps1

Write-Host "🚀 Deploying to Google Cloud Run..." -ForegroundColor Green

# Get project ID
$PROJECT_ID = gcloud config get-value project 2>$null

if (-not $PROJECT_ID) {
    Write-Host "❌ No project ID set. Run: gcloud config set project YOUR_PROJECT_ID" -ForegroundColor Red
    exit 1
}

Write-Host "📦 Project: $PROJECT_ID" -ForegroundColor Cyan
Write-Host ""

# Build and deploy
Write-Host "📦 Building container image..." -ForegroundColor Yellow
gcloud builds submit --tag "gcr.io/$PROJECT_ID/biosynth-backend" ./backend

Write-Host ""
Write-Host "🚀 Deploying to Cloud Run..." -ForegroundColor Yellow
gcloud run deploy biosynth-backend `
  --image "gcr.io/$PROJECT_ID/biosynth-backend" `
  --platform managed `
  --region us-central1 `
  --allow-unauthenticated `
  --port 8080 `
  --memory 2Gi `
  --cpu 2 `
  --timeout 300 `
  --max-instances 10 `
  --set-env-vars "PORT=8080,ENVIRONMENT=production" `
  --quiet

Write-Host ""
Write-Host "✅ Deployment complete!" -ForegroundColor Green
Write-Host ""

Write-Host "🌐 Service URL:" -ForegroundColor Cyan
$SERVICE_URL = gcloud run services describe biosynth-backend --region us-central1 --format 'value(status.url)'
Write-Host $SERVICE_URL -ForegroundColor White

Write-Host ""
Write-Host "🧪 Test health endpoint:" -ForegroundColor Cyan
Write-Host "curl $SERVICE_URL/health" -ForegroundColor White

