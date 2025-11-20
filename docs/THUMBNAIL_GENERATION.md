# Thumbnail Generation Pipeline

Complete guide for generating molecule thumbnails, 3D structures, and assets.

## Overview

The thumbnail generation system provides:
- **2D thumbnails** - Clean structural diagrams for fast loading
- **3D molfiles** - Optimized 3D coordinates for 3Dmol.js rendering
- **Base64 encoding** - Direct database storage support
- **Batch processing** - Process entire directories of molecule files

## Architecture

### Backend Services

1. **ThumbnailService** (`backend/services/thumbnail_service.py`)
   - Core service for generating thumbnails using RDKit
   - Supports both molfile and SMILES input
   - Generates 2D images and 3D structures

2. **Thumbnail API Endpoints** (`backend/routes/thumbnails.py`)
   - `/thumbnails/generate` - Returns PNG image bytes
   - `/thumbnails/generate-base64` - Returns base64-encoded image
   - `/thumbnails/generate-3d` - Generates optimized 3D molfile

3. **Batch Processing Script** (`backend/scripts/generate_molecule_assets.py`)
   - CLI tool for processing multiple molecule files
   - Generates thumbnails, 3D molfiles, and metadata

### Frontend Integration

- **API Client** (`frontend/src/lib/api.ts`)
  - `generateThumbnail()` - Get thumbnail as blob
  - `generateThumbnailBase64()` - Get thumbnail as base64
  - `generate3DMolfile()` - Get optimized 3D molfile

- **Lab Page** (`frontend/src/pages/Lab.tsx`)
  - Automatically generates thumbnails when saving molecules
  - Falls back to backend generation if canvas thumbnail unavailable

## API Usage

### Generate Thumbnail (Base64)

```typescript
import { generateThumbnailBase64 } from '../lib/api';

const result = await generateThumbnailBase64({
  smiles: 'CCO', // or molfile: '...'
  size: [600, 600] // optional
});

// result.thumbnail_b64 is a data URI: "data:image/png;base64,..."
```

### Generate 3D Molfile

```typescript
import { generate3DMolfile } from '../lib/api';

const result = await generate3DMolfile({
  smiles: 'CCO'
});

// result.molfile contains optimized 3D coordinates
```

### Direct API Calls

```bash
# Generate thumbnail (returns PNG image)
curl -X POST http://localhost:8000/thumbnails/generate \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CCO"}'

# Generate base64 thumbnail
curl -X POST http://localhost:8000/thumbnails/generate-base64 \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CCO", "size": [600, 600]}'

# Generate 3D molfile
curl -X POST http://localhost:8000/thumbnails/generate-3d \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CCO"}'
```

## Batch Processing

### Setup

1. Create input directory with molecule files:
```bash
mkdir -p molecules_raw
cp *.mol molecules_raw/
```

2. Run batch processing:
```bash
python backend/scripts/generate_molecule_assets.py \
  --input-dir molecules_raw \
  --output-dir public/molecules
```

### Output Structure

```
public/molecules/
├── aspirin.png              # 2D thumbnail
├── aspirin_3d.mol          # 3D optimized molfile
├── aspirin_metadata.json   # Metadata (SMILES, paths, etc.)
├── caffeine.png
├── caffeine_3d.mol
└── caffeine_metadata.json
```

### Options

```bash
# Skip 3D generation (faster)
python backend/scripts/generate_molecule_assets.py --no-3d

# Process specific file types
python backend/scripts/generate_molecule_assets.py \
  --extensions .mol .sdf

# Custom directories
python backend/scripts/generate_molecule_assets.py \
  --input-dir /path/to/molecules \
  --output-dir /path/to/output
```

## Supported Formats

### Input Formats
- `.mol` - MDL Molfile (V2000/V3000)
- `.sdf` - Structure Data Format
- `.pdb` - Protein Data Bank format
- `.smi` - SMILES file (one per line)

### Output Formats
- **PNG** - 2D thumbnail images (default 600x600)
- **SVG** - Vector format (optional)
- **Molfile** - 3D optimized coordinates (V2000)

## Integration with Library Page

The Library page automatically uses thumbnails:

1. **Priority order:**
   - Canvas-generated thumbnail (if available)
   - Backend-generated thumbnail (from SMILES/molfile)
   - Placeholder image

2. **Lazy loading:**
   - Thumbnails load instantly
   - 3D viewers mount only when cards enter viewport
   - Smooth fade-in transitions

3. **Auto-generation:**
   - When saving molecules in Lab, thumbnails are generated automatically
   - Missing thumbnails can be regenerated via API

## Performance Considerations

### Thumbnail Generation
- **2D thumbnails**: ~50-200ms per molecule
- **3D structures**: ~200-1000ms per molecule (depends on complexity)
- **Batch processing**: Processes sequentially to avoid memory issues

### Optimization Tips
1. **Cache thumbnails** - Store in database or CDN
2. **Pre-generate** - Use batch script for large libraries
3. **Lazy generation** - Generate on-demand when viewing
4. **Size optimization** - Use 400x400 for thumbnails, 600x600 for previews

## Error Handling

The system gracefully handles errors:

- **Invalid SMILES/molfile**: Returns 400 error with details
- **3D generation failure**: Falls back to 2D-only
- **Missing dependencies**: Clear error messages
- **Network errors**: Frontend retries with fallbacks

## Dependencies

### Backend
- `rdkit` - Chemistry toolkit (required)
- `Pillow` - Image processing (usually included with RDKit)

### Frontend
- No additional dependencies (uses existing API client)

## Troubleshooting

### RDKit Not Found
```bash
# Install RDKit (conda recommended)
conda install -c conda-forge rdkit

# Or pip (may have issues)
pip install rdkit
```

### Thumbnail Generation Fails
1. Check SMILES/molfile validity
2. Verify RDKit installation
3. Check backend logs for detailed errors

### Batch Script Issues
1. Ensure input directory exists
2. Check file permissions
3. Verify file formats are supported

## Future Enhancements

Potential improvements:
- [ ] GLB model generation (requires Blender/three.js)
- [ ] 3D snapshot rendering (headless rendering)
- [ ] Parallel batch processing
- [ ] Thumbnail caching layer
- [ ] CDN integration for asset delivery

