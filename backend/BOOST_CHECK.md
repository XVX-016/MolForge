# Boost Dependency Check

## Do You Need Boost?

### Current Setup
- Using `rdkit-pypi==2022.9.5` (pre-built wheel)
- RDKit operations: featurization, molecule processing

### Answer: **Probably NOT needed**

`rdkit-pypi` comes with **pre-compiled binaries** that include Boost internally. Most RDKit operations work fine without system Boost libraries.

---

## Test Strategy

### Step 1: Try WITHOUT Boost (Fast Build)

Use the standard `Dockerfile` (already optimized):
- ✅ Build time: **25-35 seconds**
- ✅ Image size: **~200MB**
- ✅ Only installs `gcc` and `g++` (for compiling some Python packages)

### Step 2: If RDKit Fails at Runtime

If you see errors like:
```
ImportError: libboost_system.so.1.74.0: cannot open shared object file
```

Then switch to `Dockerfile.with-boost`:
- ⏱️ Build time: **45-60 seconds**
- 📦 Image size: **~250MB**
- ✅ Includes minimal Boost runtime libraries

---

## How to Switch

### Option 1: Use Standard Dockerfile (Recommended First)
```bash
# Already set up - just deploy
```

### Option 2: Use Boost Version
In Railway Settings → Build:
- **Dockerfile Path**: `Dockerfile.with-boost`

---

## What Boost Libraries Are Actually Needed?

If you need Boost, you only need **runtime libraries**, not dev packages:

**Minimal (what we use):**
- `libboost-system1.74.0` - System utilities
- `libboost-thread1.74.0` - Threading support

**NOT needed:**
- `libboost-dev` - Full dev package (500MB+)
- `build-essential` - Only needed for compilation
- Other Boost modules - RDKit doesn't use them

---

## Build Time Comparison

| Version | Build Time | Image Size | Boost |
|---------|-----------|------------|-------|
| **Standard** (no Boost) | 25-35s | ~200MB | ❌ |
| **With Boost** (minimal) | 45-60s | ~250MB | ✅ Runtime only |
| **Full Boost Dev** (old) | 2-3 min | ~700MB | ✅ Full dev |

---

## Recommendation

1. **Start with standard Dockerfile** (no Boost)
2. **Deploy and test** RDKit operations
3. **If errors occur**, switch to `Dockerfile.with-boost`
4. **Monitor logs** for Boost-related errors

---

## RDKit Usage in Your Backend

Your backend uses RDKit in:
- `backend/ai/featurizer.py` - SMILES featurization
- `backend/ai/train_property_predictor.py` - Training data processing
- `backend/simulation/reaction_engine.py` - Reaction processing

All of these should work with `rdkit-pypi` without system Boost.

---

## Quick Test

After deployment, test RDKit import:
```python
# In Railway logs or test endpoint
import rdkit
from rdkit import Chem
print("RDKit OK!")
```

If this works, you don't need Boost! ✅

