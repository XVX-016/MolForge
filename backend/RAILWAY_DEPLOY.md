# Railway Deployment Guide

## Quick Start

### 1. Create Railway Account

1. Go to [railway.app](https://railway.app)
2. Sign up with GitHub
3. Create a new project

### 2. Deploy from GitHub

**Option A: Auto-deploy from GitHub (Recommended)**

1. In Railway dashboard, click **"New Project"**
2. Select **"Deploy from GitHub repo"**
3. Choose your `biosynth-monorepo` repository
4. Railway will auto-detect the `backend/Dockerfile`
5. Set the **Root Directory** to `backend`:
   - Go to Settings → Source
   - Set Root Directory: `backend`

**Option B: Deploy from CLI**

```bash
# Install Railway CLI
npm i -g @railway/cli

# Login
railway login

# Initialize project
cd backend
railway init

# Deploy
railway up
```

### 3. Configure Environment Variables

In Railway dashboard → Your Service → Variables, add:

```env
# Database (Railway provides PostgreSQL automatically)
# Or use SQLite for development
DATABASE_URL=sqlite:///./biosynth.db

# CORS - Add your frontend URLs (comma-separated)
CORS_ORIGINS=https://your-app.web.app,https://your-app.firebaseapp.com,http://localhost:5173

# Optional: Model paths (if you upload model files)
MODEL_WEIGHTS_PATH=backend/weights/property_predictor.pt
ONNX_MODEL_PATH=backend/weights/predictor.onnx

# Optional: Security
SECRET_KEY=your-secret-key-here
ENVIRONMENT=production
```

### 4. Get Your Railway URL

After deployment:
1. Go to your service in Railway
2. Click on the service
3. Copy the **Public Domain** (e.g., `biosynth-backend-production.up.railway.app`)

### 5. Update Frontend API URL

Update your frontend to use the Railway URL:

**Option A: Environment Variable (Recommended)**

Create `.env` in `frontend/`:
```env
VITE_API_URL=https://your-railway-url.up.railway.app
```

**Option B: Update `frontend/src/lib/api.ts`**

```typescript
const API_URL = import.meta.env.VITE_API_URL || 
  'https://your-railway-url.up.railway.app'
```

### 6. Update CORS in Railway

In Railway → Variables, set:
```
CORS_ORIGINS=https://your-firebase-app.web.app,https://your-firebase-app.firebaseapp.com,http://localhost:5173
```

---

## Railway-Specific Configuration

### Port Configuration

Railway automatically provides a `PORT` environment variable. The Dockerfile is configured to use it:

```dockerfile
CMD uvicorn app:app --host 0.0.0.0 --port ${PORT:-8000}
```

### Database Options

**Option 1: Railway PostgreSQL (Recommended for Production)**

1. In Railway, click **"New"** → **"Database"** → **"Add PostgreSQL"**
2. Railway will automatically set `DATABASE_URL` environment variable
3. Update your backend to use PostgreSQL connection string

**Option 2: SQLite (Development/Testing)**

Keep using SQLite by setting:
```
DATABASE_URL=sqlite:///./biosynth.db
```

Note: SQLite files are ephemeral on Railway and will be lost on redeploy.

### Build Configuration

Railway will automatically:
- Detect `backend/Dockerfile`
- Build the Docker image
- Deploy it

The `railway.json` file provides additional configuration if needed.

---

## Testing Your Deployment

### 1. Test Health Endpoint

```bash
curl https://your-railway-url.up.railway.app/health
# Should return: {"status":"healthy"}
```

### 2. Test API Root

```bash
curl https://your-railway-url.up.railway.app/
# Should return: {"message":"MolForge Backend API","version":"0.1.0"}
```

### 3. Test from Frontend

Make sure your frontend can connect:
- Check browser console for CORS errors
- Test API calls from your React app

---

## Troubleshooting

### Issue: Build Fails

**Solution:**
- Check Railway logs for build errors
- Ensure `requirements.txt` is in `backend/` directory
- Verify Dockerfile is correct

### Issue: CORS Errors

**Solution:**
- Add your frontend URL to `CORS_ORIGINS` in Railway variables
- Format: `https://your-app.web.app,http://localhost:5173` (comma-separated)

### Issue: Database Connection Fails

**Solution:**
- If using PostgreSQL, ensure Railway PostgreSQL service is running
- Check `DATABASE_URL` is set correctly in Railway variables
- For SQLite, ensure path is correct

### Issue: Port Already in Use

**Solution:**
- Railway handles ports automatically
- Don't hardcode port 8000 - use `$PORT` environment variable
- The Dockerfile is already configured for this

### Issue: Module Not Found

**Solution:**
- Check `requirements.txt` includes all dependencies
- Verify Python version in Dockerfile matches your local (3.11)

---

## Railway Pricing

- **Free Tier**: $5 credit/month
- **Hobby Plan**: $5/month (includes $5 credit)
- **Pro Plan**: $20/month

For most projects, the free tier is sufficient for development.

---

## Continuous Deployment

Railway automatically deploys when you push to your connected branch (usually `main` or `master`).

To disable auto-deploy:
1. Go to Settings → Source
2. Toggle off "Auto Deploy"

---

## Custom Domain (Optional)

1. Go to Settings → Networking
2. Click "Generate Domain" or add custom domain
3. Update `CORS_ORIGINS` to include your custom domain

---

## Monitoring

Railway provides:
- Real-time logs
- Metrics (CPU, Memory, Network)
- Deployment history

Access via Railway dashboard → Your Service → Metrics

---

## Next Steps

1. ✅ Deploy backend to Railway
2. ✅ Get Railway URL
3. ✅ Update frontend API URL
4. ✅ Update CORS settings
5. ✅ Test API endpoints
6. ✅ Deploy frontend to Firebase Hosting
7. ✅ Update Firebase frontend to use Railway backend URL

---

## Example Railway Variables

```
DATABASE_URL=postgresql://postgres:password@containers-us-west-xxx.railway.app:5432/railway
CORS_ORIGINS=https://biosynth-app.web.app,https://biosynth-app.firebaseapp.com,http://localhost:5173
ENVIRONMENT=production
SECRET_KEY=your-secret-key-here
LOG_LEVEL=INFO
```

---

## Support

- Railway Docs: https://docs.railway.app
- Railway Discord: https://discord.gg/railway
- Railway Status: https://status.railway.app

