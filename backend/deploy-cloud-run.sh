#!/bin/bash
# Quick deploy script for Google Cloud Run
# Usage: ./deploy-cloud-run.sh

set -e

# Get project ID
PROJECT_ID=$(gcloud config get-value project 2>/dev/null)

if [ -z "$PROJECT_ID" ]; then
    echo "❌ No project ID set. Run: gcloud config set project YOUR_PROJECT_ID"
    exit 1
fi

echo "🚀 Deploying to Google Cloud Run..."
echo "📦 Project: $PROJECT_ID"
echo ""

# Build and deploy
echo "📦 Building container image..."
gcloud builds submit --tag gcr.io/$PROJECT_ID/biosynth-backend ./backend

echo ""
echo "🚀 Deploying to Cloud Run..."
gcloud run deploy biosynth-backend \
  --image gcr.io/$PROJECT_ID/biosynth-backend \
  --platform managed \
  --region us-central1 \
  --allow-unauthenticated \
  --port 8080 \
  --memory 2Gi \
  --cpu 2 \
  --timeout 300 \
  --max-instances 10 \
  --set-env-vars "PORT=8080,ENVIRONMENT=production" \
  --quiet

echo ""
echo "✅ Deployment complete!"
echo ""
echo "🌐 Service URL:"
gcloud run services describe biosynth-backend \
  --region us-central1 \
  --format 'value(status.url)'

echo ""
echo "🧪 Test health endpoint:"
SERVICE_URL=$(gcloud run services describe biosynth-backend --region us-central1 --format 'value(status.url)')
echo "curl $SERVICE_URL/health"

