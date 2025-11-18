# Railway Deployment - Quick Start Checklist

## ✅ Pre-Deployment Checklist

- [x] Dockerfile updated for Railway (uses `$PORT` env var)
- [x] CORS configuration supports environment variables
- [x] Railway configuration files created
- [ ] Code pushed to GitHub
- [ ] Railway account created

---

## 🚀 Step-by-Step Deployment

### 1. Push Code to GitHub

```bash
git add .
git commit -m "Prepare for Railway deployment"
git push origin main
```

### 2. Create Railway Project

1. Go to [railway.app](https://railway.app)
2. Click **"New Project"**
3. Select **"Deploy from GitHub repo"**
4. Choose your `biosynth-monorepo` repository
5. Railway will detect the Dockerfile

### 3. Configure Root Directory

1. In Railway dashboard → Your Service → Settings
2. Go to **"Source"** section
3. Set **Root Directory**: `backend`
4. Save

### 4. Set Environment Variables

Go to Railway → Your Service → Variables, add:

```env
CORS_ORIGINS=https://your-firebase-app.web.app,https://your-firebase-app.firebaseapp.com,http://localhost:5173
ENVIRONMENT=production
```

**Note:** Replace `your-firebase-app` with your actual Firebase hosting URL.

### 5. Deploy

Railway will automatically:
- Build the Docker image
- Deploy it
- Provide a public URL

### 6. Get Your Backend URL

1. Railway dashboard → Your Service
2. Click on the service
3. Copy the **Public Domain** (e.g., `biosynth-backend-production.up.railway.app`)

### 7. Update Frontend

**Option A: Environment Variable (Recommended)**

Create `frontend/.env.production`:
```env
VITE_API_URL=https://your-railway-url.up.railway.app
```

**Option B: Update `frontend/src/lib/api.ts`**

```typescript
const API_URL = import.meta.env.VITE_API_URL || 
  'https://your-railway-url.up.railway.app'
```

### 8. Test Deployment

```bash
# Test health endpoint
curl https://your-railway-url.up.railway.app/health

# Should return: {"status":"healthy"}
```

### 9. Deploy Frontend to Firebase

```bash
cd frontend
npm run build
firebase deploy --only hosting
```

---

## 🔧 Troubleshooting

### Build Fails
- Check Railway logs
- Ensure `requirements.txt` is in `backend/` directory
- Verify Python 3.11 in Dockerfile

### CORS Errors
- Add frontend URL to `CORS_ORIGINS` in Railway variables
- Format: `https://app.web.app,http://localhost:5173` (comma-separated, no spaces)

### 500 Errors
- Check Railway logs for Python errors
- Verify all dependencies in `requirements.txt`
- Check database connection if using PostgreSQL

---

## 📝 Files Created/Modified

- ✅ `backend/Dockerfile` - Updated for Railway PORT
- ✅ `backend/config.py` - CORS from environment variables
- ✅ `backend/railway.json` - Railway configuration
- ✅ `backend/railway.toml` - Alternative Railway config
- ✅ `backend/.railwayignore` - Files to exclude
- ✅ `backend/RAILWAY_DEPLOY.md` - Full deployment guide

---

## 🎯 Next Steps After Deployment

1. ✅ Test all API endpoints
2. ✅ Update frontend to use Railway URL
3. ✅ Deploy frontend to Firebase
4. ✅ Test full stack integration
5. ✅ Monitor Railway logs for errors

---

## 💰 Railway Pricing

- **Free Tier**: $5 credit/month (usually enough for development)
- **Hobby**: $5/month
- **Pro**: $20/month

Most projects can start on the free tier!

---

## 📚 Resources

- Full Guide: `backend/RAILWAY_DEPLOY.md`
- Railway Docs: https://docs.railway.app
- Railway Status: https://status.railway.app

