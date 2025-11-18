# Railway Configuration Verification

## ✅ Critical Settings Checklist

### 1. Root Directory
**Railway Dashboard → Settings → Deployment**
- **Root Directory**: `backend` (NOT `/backend`, NOT `./backend`)
- Must be exactly: `backend`

### 2. Dockerfile Path
**Railway Dashboard → Settings → Build**
- **Dockerfile Path**: `Dockerfile` (relative to root directory)
- OR: Leave as auto-detect

### 3. Start Command
**Railway Dashboard → Settings → Deploy**
- **Start Command**: `uvicorn app:app --host 0.0.0.0 --port $PORT`
- OR: Leave empty (uses Dockerfile CMD)

### 4. Builder Type
**Railway Dashboard → Settings → Build**
- **Builder**: `Dockerfile` (NOT Nix, NOT Buildpacks)
- Turn OFF auto-detection if it's using wrong builder

---

## 🔍 How to Verify Railway is Using Correct Dockerfile

### Check Build Logs

After deployment, Railway build logs MUST show:

```
=== Using Detected Dockerfile backend/Dockerfile ===
FROM rdkit/rdkit:latest-py3.11
```

**If you see:**
```
FROM python:3.11-slim
RUN apt-get install build-essential
```

**Then Railway is NOT using your Dockerfile!**

---

## 🚨 Common Issues

### Issue 1: Dockerfile Not Committed to Git

**Check:**
```bash
git status
git log --oneline --all -- backend/Dockerfile
```

**Fix:**
```bash
git add backend/Dockerfile backend/requirements-model.txt
git commit -m "Add optimized Dockerfile with RDKit base image"
git push
```

### Issue 2: Wrong Root Directory

**Wrong:**
- `/backend`
- `./backend`
- `backend/`
- `/workspace/backend`

**Correct:**
- `backend` (exactly this)

### Issue 3: Wrong Dockerfile Path

**Wrong:**
- `/backend/Dockerfile`
- `./Dockerfile`
- `/app/Dockerfile`

**Correct:**
- `Dockerfile` (relative to root directory)

---

## ✅ Verification Steps

1. **Check Git Status**
   ```bash
   git status
   # Should show Dockerfile as tracked
   ```

2. **Verify Dockerfile Content**
   ```bash
   head -5 backend/Dockerfile
   # Should show: FROM rdkit/rdkit:latest-py3.11
   ```

3. **Check Railway Build Logs**
   - First line should show: `FROM rdkit/rdkit:latest-py3.11`
   - NOT: `FROM python:3.11-slim`

4. **Check Image Size**
   - Should be ~1GB (not 6.4GB)

---

## 🎯 Expected Build Log Output

**Correct (using RDKit image):**
```
Step 1/10 : FROM rdkit/rdkit:latest-py3.11
Step 2/10 : WORKDIR /app
Step 3/10 : COPY requirements-model.txt .
Step 4/10 : RUN pip install torch==2.4.1+cpu ...
```

**Wrong (old Dockerfile):**
```
Step 1/15 : FROM python:3.11-slim
Step 2/15 : RUN apt-get update
Step 3/15 : RUN apt-get install build-essential libboost-dev
```

---

## 📝 Files Created

- ✅ `backend/Dockerfile` - Optimized with RDKit base image
- ✅ `backend/requirements-model.txt` - Heavy ML dependencies (cached separately)
- ✅ `backend/requirements.txt` - Lightweight dependencies
- ✅ `backend/.dockerignore` - Excludes unnecessary files

---

## 🚀 Next Steps

1. ✅ Verify Dockerfile is committed to Git
2. ✅ Set Railway Root Directory to `backend`
3. ✅ Verify Dockerfile Path is `Dockerfile`
4. ✅ Push changes
5. ✅ Check build logs for `FROM rdkit/rdkit:latest-py3.11`
6. ✅ Verify image size is ~1GB (not 6.4GB)

---

**If build logs still show `python:3.11-slim`, Railway is NOT using your Dockerfile!**

