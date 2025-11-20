# Firebase Library Scaffold - Summary

This document summarizes all the files created for the Firebase + Library integration.

## ✅ Files Created

### Core Firebase Files

1. **`src/firebase.ts`**
   - Firebase initialization with Auth + Firestore
   - Exports `auth` and `db` for use throughout the app

2. **`src/lib/firebaseMoleculeStore.ts`**
   - Firestore CRUD helpers:
     - `saveMolecule(userId, molecule)` - Save a molecule
     - `listMolecules(userId)` - List all molecules for a user
     - `deleteMolecule(userId, moleculeId)` - Delete a molecule
     - `getMolecule(userId, moleculeId)` - Get a single molecule
     - `searchMolecules(userId, query)` - Search molecules by name/SMILES

3. **`src/lib/smiles.ts`**
   - SMILES utility functions:
     - `isValidSMILES()` - Basic SMILES validation
     - `normalizeSMILES()` - Normalize SMILES strings
     - `searchBySMILES()` - Search molecules by SMILES pattern

### UI Components

4. **`src/components/ui/GlassCard.tsx`**
   - Glassmorphic card component with backdrop blur
   - Reusable for any glass-style UI element

5. **`src/components/molecules/MiniMoleculeViewer.tsx`**
   - Mini 3D molecule viewer component
   - Can be embedded in molecule cards
   - Uses React Three Fiber for rendering

### Pages

6. **`src/pages/LibraryPage.tsx`**
   - Alternative Library page using Firebase instead of backend API
   - Includes authentication check
   - Full search and pagination support

### Styles

7. **`src/styles/glass.css`**
   - Glassmorphic UI styles
   - Includes `.glass-card`, `.glass-button`, `.glass-input` classes
   - Automatically imported via `index.css`

### Configuration Files

8. **`firebase.rules`**
   - Firestore security rules
   - Users can only access their own molecules
   - Deploy with: `firebase deploy --only firestore:rules`

9. **`README_FIREBASE.md`**
   - Complete setup guide
   - Step-by-step Firebase configuration
   - Troubleshooting tips

## 📦 Dependencies

Firebase is already installed in `package.json`:
```json
"firebase": "^12.6.0"
```

## 🚀 Quick Start

### 1. Set up Firebase

Follow the instructions in `README_FIREBASE.md` to:
- Create a Firebase project
- Enable Firestore
- Get your configuration credentials

### 2. Configure Environment Variables

Create `frontend/.env` (copy from `.env.example`):
```env
VITE_FIREBASE_API_KEY=your-api-key
VITE_FIREBASE_AUTH_DOMAIN=your-project.firebaseapp.com
VITE_FIREBASE_PROJECT_ID=your-project-id
VITE_FIREBASE_STORAGE_BUCKET=your-project.appspot.com
VITE_FIREBASE_MESSAGING_SENDER_ID=your-sender-id
VITE_FIREBASE_APP_ID=your-app-id
```

### 3. Deploy Firestore Rules

```bash
cd frontend
firebase deploy --only firestore:rules
```

### 4. Use in Your Code

#### Save a molecule:
```typescript
import { auth } from '../firebase';
import { saveMolecule } from '../lib/firebaseMoleculeStore';

const userId = auth.currentUser?.uid;
if (userId) {
  await saveMolecule(userId, {
    name: 'Benzene',
    smiles: 'c1ccccc1',
    formula: 'C6H6',
  });
}
```

#### List molecules:
```typescript
import { listMolecules } from '../lib/firebaseMoleculeStore';

const molecules = await listMolecules(userId);
```

#### Use glass components:
```tsx
import GlassCard from '../components/ui/GlassCard';

<GlassCard className="p-6">
  <h2>My Content</h2>
</GlassCard>
```

## 🔄 Integration Options

### Option 1: Use LibraryPage (Firebase-only)

Update your routing to use `LibraryPage` instead of `Library`:
```tsx
import LibraryPage from './pages/LibraryPage';
// Use LibraryPage in your routes
```

### Option 2: Hybrid Approach

Keep using `Library.tsx` (backend API) but add Firebase as an option:
- Add a toggle to switch between backend and Firebase
- Or use Firebase for authenticated users, backend for guests

### Option 3: Update Existing Library.tsx

Modify `Library.tsx` to use Firebase:
```typescript
// Replace: import { listMolecules } from '../lib/api';
import { listMolecules } from '../lib/firebaseMoleculeStore';
import { auth } from '../firebase';

// Then use:
const userId = auth.currentUser?.uid;
if (userId) {
  const molecules = await listMolecules(userId);
}
```

## 🎨 Glassmorphic UI

The glass styles are automatically available. Use them like:

```tsx
<div className="glass-card">
  Frosted glass effect
</div>

<button className="glass-button">
  Glass button
</button>

<input className="glass-input" />
```

## 📝 Next Steps

1. **Add Authentication UI**
   - Create login/signup components
   - Wire up Firebase Auth

2. **Update Lab Page**
   - Add "Save to Firebase" button
   - Use `saveMolecule()` from `firebaseMoleculeStore`

3. **Enhance MoleculeCard**
   - Add `MiniMoleculeViewer` for 3D previews
   - Apply glass styles

4. **Add User Profile**
   - Show user's molecule count
   - Add export/import functionality

## 🔍 File Structure

```
frontend/
├── src/
│   ├── firebase.ts                    # Firebase init
│   ├── lib/
│   │   ├── firebaseMoleculeStore.ts   # Firestore CRUD
│   │   └── smiles.ts                 # SMILES utilities
│   ├── components/
│   │   ├── ui/
│   │   │   └── GlassCard.tsx         # Glass card component
│   │   └── molecules/
│   │       └── MiniMoleculeViewer.tsx # Mini 3D viewer
│   ├── pages/
│   │   └── LibraryPage.tsx           # Firebase-based Library
│   └── styles/
│       └── glass.css                 # Glassmorphic styles
├── firebase.rules                     # Firestore security rules
├── README_FIREBASE.md                 # Setup guide
└── .env.example                      # Environment template
```

## 🐛 Troubleshooting

See `README_FIREBASE.md` for detailed troubleshooting.

Common issues:
- **Missing env vars**: Make sure `.env` exists and has all Firebase config
- **Permission errors**: Deploy Firestore rules
- **Auth errors**: Enable Authentication in Firebase Console

## 📚 Additional Resources

- [Firebase Documentation](https://firebase.google.com/docs)
- [Firestore Documentation](https://firebase.google.com/docs/firestore)
- [Firebase Auth Documentation](https://firebase.google.com/docs/auth)

