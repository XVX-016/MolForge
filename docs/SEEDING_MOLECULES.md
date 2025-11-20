# Seeding Molecules for Library

Guide for generating and seeding sample molecules into your library.

## Overview

Two methods are available to generate sample molecules:

1. **Backend Python Script** - Generates molecules with all properties (molfile, thumbnails)
2. **Frontend UI** - Interactive page to seed molecules via the web interface

## Method 1: Backend Python Script

### Usage

```bash
# Generate sample molecules JSON file
python backend/scripts/generate_sample_molecules.py

# Generate with SQL output
python backend/scripts/generate_sample_molecules.py --sql --user-id YOUR_USER_ID

# Custom output file
python backend/scripts/generate_sample_molecules.py --output my_molecules.json
```

### What It Generates

For each molecule, the script creates:
- ✅ **SMILES** string
- ✅ **Formula** (e.g., C9H8O4)
- ✅ **3D Molfile** (optimized coordinates)
- ✅ **2D Thumbnail** (base64 PNG image)
- ✅ **Properties** (stability, toxicity, solubility, etc.)

### Output

1. **JSON File** (`sample_molecules.json`)
   - Contains all molecule data
   - Can be imported programmatically
   - Useful for testing and development

2. **SQL Statements** (with `--sql` flag)
   - Ready-to-run INSERT statements
   - Can be pasted into Supabase SQL editor
   - Replace `YOUR_USER_ID` with actual user ID

### Sample Molecules Included

- Aspirin (C9H8O4)
- Caffeine (C8H10N4O2)
- Glucose (C6H12O6)
- Ethanol (C2H6O)
- Paracetamol (C8H9NO2)
- Ibuprofen (C13H18O2)
- Vitamin C (C6H8O6)
- Nicotine (C10H14N2)
- Serotonin (C10H12N2O)
- Dopamine (C8H11NO2)
- And 10 more...

## Method 2: Frontend UI

### Access

Navigate to: `/seed-library` in your application

### Features

- ✅ Interactive UI with progress tracking
- ✅ Authentication check (requires sign-in)
- ✅ Real-time progress bar
- ✅ Success/failure reporting
- ✅ Automatic generation of molfiles and thumbnails

### Usage Steps

1. **Sign in** to your account
2. Navigate to `/seed-library`
3. Review the list of sample molecules
4. Click **"Seed Library"** button
5. Wait for generation (may take a few minutes)
6. Visit `/library` to see your new molecules

### What Happens

1. For each molecule:
   - Calls backend API to generate 3D molfile
   - Calls backend API to generate thumbnail
   - Saves to Supabase with all properties
2. Progress is shown in real-time
3. Results summary displayed at completion

## Sample Molecule Data

Each generated molecule includes:

```typescript
{
  name: "Aspirin",
  smiles: "CC(=O)OC1=CC=CC=C1C(=O)O",
  formula: "C9H8O4",
  molfile: "...", // 3D optimized molfile
  thumbnail_b64: "data:image/png;base64,...", // 2D structure image
  properties: {
    stability: 0.85,
    toxicity: 0.15,
    solubility: 0.72,
    bioavailability: 0.78,
    novelty: 0.45
  }
}
```

## Requirements

### Backend Script
- Python 3.8+
- RDKit installed (`pip install rdkit` or `conda install rdkit`)
- Backend services available (thumbnail_service)

### Frontend UI
- Backend API running
- Supabase configured
- User authentication set up

## Troubleshooting

### Backend Script Issues

**RDKit not found:**
```bash
# Install via conda (recommended)
conda install -c conda-forge rdkit

# Or via pip (may have issues)
pip install rdkit-pypi
```

**Generation fails:**
- Check that backend services are accessible
- Verify RDKit installation
- Check molecule SMILES validity

### Frontend UI Issues

**"Authentication Required":**
- Sign in via `/supabase-test` or your auth page
- Ensure Supabase is configured

**Generation fails:**
- Check backend API is running
- Verify API endpoints are accessible
- Check browser console for errors

**Slow generation:**
- Normal - each molecule requires API calls
- 20 molecules ≈ 2-5 minutes
- Progress bar shows real-time status

## Customization

### Add More Molecules

**Backend Script:**
Edit `backend/scripts/generate_sample_molecules.py`:
```python
SAMPLE_MOLECULES.append({
    "name": "Your Molecule",
    "smiles": "YOUR_SMILES",
    "formula": "CxHyOz",
    "description": "Description"
})
```

**Frontend:**
Edit `frontend/src/utils/seedMolecules.ts`:
```typescript
export const SAMPLE_MOLECULES: SampleMolecule[] = [
  // Add your molecules here
  {
    name: 'Your Molecule',
    smiles: 'YOUR_SMILES',
    formula: 'CxHyOz',
    description: 'Description',
  },
];
```

## Best Practices

1. **Test First**: Seed a few molecules before seeding all
2. **Backup**: Export existing molecules before seeding
3. **User-Specific**: Each user should seed their own library
4. **Cleanup**: Delete test molecules after development

## Next Steps

After seeding:
1. Visit `/library` to see your molecules
2. Test search functionality
3. Verify 3D viewers load correctly
4. Check thumbnails display properly

## API Endpoints Used

- `POST /convert/smiles` - Generate 3D molfile
- `POST /thumbnails/generate-base64` - Generate thumbnail
- Supabase `molecules` table - Store molecule data

