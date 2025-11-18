# Google Cloud Run - Quick Start

## 🚀 5-Minute Deployment

### Step 1: Install Google Cloud SDK

**Windows:**
```powershell
# Download installer
(New-Object Net.WebClient).DownloadFile("https://dl.google.com/dl/cloudsdk/channels/rapid/GoogleCloudSDKInstaller.exe", "$env:Temp\GoogleCloudSDKInstaller.exe")
& $env:Temp\GoogleCloudSDKInstaller.exe
```

Or download from: https://cloud.google.com/sdk/docs/install

### Step 2: Authenticate

```bash
gcloud auth login
gcloud auth application-default login
```

### Step 3: Create/Select Project

```bash
# Create new project
gcloud projects create biosynth-backend --name="BioSynth Backend"

# Set as current project
gcloud config set project biosynth-backend
```

### Step 4: Enable APIs

```bash
gcloud services enable cloudbuild.googleapis.com
gcloud services enable run.googleapis.com
gcloud services enable containerregistry.googleapis.com
```

### Step 5: Deploy

**Option A: Using PowerShell Script (Windows)**
```powershell
cd backend
.\deploy-cloud-run.ps1
```

**Option B: Manual Deploy**
```bash
# From repo root
export PROJECT_ID=$(gcloud config get-value project)
gcloud builds submit --tag gcr.io/$PROJECT_ID/biosynth-backend ./backend

gcloud run deploy biosynth-backend \
  --image gcr.io/$PROJECT_ID/biosynth-backend \
  --platform managed \
  --region us-central1 \
  --allow-unauthenticated \
  --port 8080 \
  --memory 2Gi \
  --cpu 2 \
  --set-env-vars "PORT=8080,ENVIRONMENT=production"
```

### Step 6: Get Your URL

```bash
gcloud run services describe biosynth-backend \
  --region us-central1 \
  --format 'value(status.url)'
```

Copy this URL and update your frontend!

---

## ✅ Configuration Checklist

- [ ] Google Cloud SDK installed
- [ ] Authenticated (`gcloud auth login`)
- [ ] Project created/selected
- [ ] APIs enabled
- [ ] Backend deployed
- [ ] Service URL copied
- [ ] Frontend API URL updated
- [ ] CORS_ORIGINS set

---

## 🔧 Set Environment Variables

```bash
gcloud run services update biosynth-backend \
  --region us-central1 \
  --update-env-vars "CORS_ORIGINS=https://your-app.web.app,http://localhost:5173"
```

---

## 📝 Update Frontend

In `frontend/src/lib/api.ts`:

```typescript
const API_URL = import.meta.env.VITE_API_URL || 
  'https://biosynth-backend-xxxxx-uc.a.run.app'
```

---

## 🧪 Test Deployment

```bash
# Get service URL
SERVICE_URL=$(gcloud run services describe biosynth-backend --region us-central1 --format 'value(status.url)')

# Test health
curl $SERVICE_URL/health
# Should return: {"status":"healthy"}
```

---

## 💰 Cost

**Free Tier:**
- 2 million requests/month
- 360,000 GB-seconds
- 180,000 vCPU-seconds

**After Free Tier:**
- ~$0.40 per million requests
- Very affordable for most apps

---

## 📚 Full Documentation

See `backend/CLOUD_RUN_DEPLOY.md` for complete guide.

---

**Your backend is ready for Cloud Run! 🚀**

