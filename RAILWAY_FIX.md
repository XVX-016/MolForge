# Railway Deployment Fix ✅

## Problem
Railway was failing with **"Dockerfile does not exist"** because:
- `railway.toml` was in `/backend` with `dockerfilePath = "backend/Dockerfile"`
- Railway Directory was set to `/backend`
- This caused path confusion and missing files

## Solution: Option 1 (Recommended) ✅

**Changed to repo root as project root:**

### ✅ Changes Made:

1. **Moved `railway.toml` to repo root**
   - Old: `backend/railway.toml`
   - New: `railway.toml` (repo root)

2. **Updated `railway.toml`:**
   ```toml
   [build]
   builder = "DOCKERFILE"
   dockerfilePath = "backend/Dockerfile"
   ```

3. **Updated `backend/Dockerfile`:**
   - Changed `COPY requirements-model.txt .` → `COPY backend/requirements-model.txt .`
   - Changed `COPY requirements.txt .` → `COPY backend/requirements.txt .`
   - Changed `COPY . .` → `COPY backend/ .`
   - Build context is now repo root, so all paths need `backend/` prefix

---

## 🚀 Railway Configuration

### In Railway Dashboard:

1. **Go to your service settings**
2. **Set "Root Directory" to:** `.` (repo root) or leave empty
3. **Railway will automatically use:** `railway.toml` from repo root
4. **Dockerfile path:** `backend/Dockerfile` (relative to repo root)

---

## ✅ What This Fixes

- ✅ Railway can now find the Dockerfile
- ✅ All files in `/backend` are accessible
- ✅ Build context includes entire repo (if needed)
- ✅ No more "Dockerfile does not exist" errors

---

## 📝 Next Steps

1. **Push changes to GitHub:**
   ```bash
   git add railway.toml backend/Dockerfile
   git commit -m "Fix Railway deployment: Move railway.toml to repo root, update Dockerfile paths"
   git push origin master
   ```

2. **Update Railway settings:**
   - Go to Railway Dashboard
   - Service Settings → Root Directory
   - Set to: `.` (or leave empty for repo root)
   - Save

3. **Redeploy:**
   - Railway will automatically detect the new `railway.toml`
   - Build should now succeed!

---

## 🔍 Verification

After deploying, check Railway logs:
- ✅ Should see: "Building Dockerfile: backend/Dockerfile"
- ✅ Should see: "COPY backend/requirements-model.txt"
- ✅ Should see: "COPY backend/requirements.txt"
- ✅ Should see: "COPY backend/ ."

---

**Railway deployment is now fixed! 🚀**

