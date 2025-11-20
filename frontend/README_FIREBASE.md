# Firebase Setup Guide

This guide will help you set up Firebase for the Biosynth Molecule Library.

## Prerequisites

- A Firebase account (sign up at [firebase.google.com](https://firebase.google.com))
- Node.js and npm installed

## Step 1: Create a Firebase Project

1. Go to [Firebase Console](https://console.firebase.google.com)
2. Click "Add project" or select an existing project
3. Follow the setup wizard:
   - Enter a project name (e.g., "biosynth-library")
   - Enable Google Analytics (optional)
   - Click "Create project"

## Step 2: Enable Firestore Database

1. In your Firebase project, go to **Build** → **Firestore Database**
2. Click "Create database"
3. Choose **Start in test mode** (we'll secure it with rules later)
4. Select a location for your database
5. Click "Enable"

## Step 3: Enable Firebase Authentication (Optional but Recommended)

1. Go to **Build** → **Authentication**
2. Click "Get started"
3. Enable **Email/Password** sign-in method
4. Save

## Step 4: Get Your Firebase Configuration

1. Go to **Project Settings** (gear icon) → **General**
2. Scroll down to "Your apps"
3. Click the **Web** icon (`</>`) to add a web app
4. Register your app with a nickname (e.g., "Biosynth Frontend")
5. Copy the Firebase configuration object

## Step 5: Configure Environment Variables

1. Copy `.env.example` to `.env` in the `frontend/` directory:
   ```bash
   cp frontend/.env.example frontend/.env
   ```

2. Open `frontend/.env` and fill in your Firebase credentials:
   ```env
   VITE_FIREBASE_API_KEY=your-api-key-here
   VITE_FIREBASE_AUTH_DOMAIN=your-project-id.firebaseapp.com
   VITE_FIREBASE_PROJECT_ID=your-project-id
   VITE_FIREBASE_STORAGE_BUCKET=your-project-id.appspot.com
   VITE_FIREBASE_MESSAGING_SENDER_ID=your-messaging-sender-id
   VITE_FIREBASE_APP_ID=your-app-id
   ```

## Step 6: Deploy Firestore Security Rules

1. Install Firebase CLI (if not already installed):
   ```bash
   npm install -g firebase-tools
   ```

2. Login to Firebase:
   ```bash
   firebase login
   ```

3. Initialize Firebase in your project (if not already done):
   ```bash
   cd frontend
   firebase init firestore
   ```
   - Select your Firebase project
   - Use `firebase.rules` as your rules file

4. Deploy the security rules:
   ```bash
   firebase deploy --only firestore:rules
   ```

## Step 7: Install Dependencies

Firebase is already included in `package.json`. If you need to reinstall:

```bash
cd frontend
npm install
```

## Step 8: Test the Setup

1. Start your development server:
   ```bash
   npm run dev
   ```

2. Navigate to the Library page
3. Try saving a molecule (you'll need to implement authentication first)

## Data Structure

Molecules are stored in Firestore with the following structure:

```
users/
  {userId}/
    molecules/
      {moleculeId}/
        name: string
        smiles?: string
        formula?: string
        json_graph?: string
        properties?: string
        thumbnail_b64?: string
        userId: string
        createdAt: number (timestamp)
```

## Security Rules

The `firebase.rules` file enforces:
- Users can only read/write their own molecules
- All other access is denied

**Important**: Make sure to deploy these rules before going to production!

## Using Firebase in Your Code

### Import Firebase services:

```typescript
import { auth, db } from '../firebase';
import { saveMolecule, listMolecules, deleteMolecule } from '../lib/firebaseMoleculeStore';
```

### Save a molecule:

```typescript
const userId = auth.currentUser?.uid;
if (userId) {
  const moleculeId = await saveMolecule(userId, {
    name: 'Benzene',
    smiles: 'c1ccccc1',
    formula: 'C6H6',
  });
}
```

### List molecules:

```typescript
const userId = auth.currentUser?.uid;
if (userId) {
  const molecules = await listMolecules(userId);
}
```

### Delete a molecule:

```typescript
const userId = auth.currentUser?.uid;
if (userId && moleculeId) {
  await deleteMolecule(userId, moleculeId);
}
```

## Troubleshooting

### "Firebase: Error (auth/configuration-not-found)"
- Make sure your `.env` file exists and has all required variables
- Restart your dev server after changing `.env`

### "Missing or insufficient permissions"
- Check that your Firestore rules are deployed
- Verify that the user is authenticated
- Check that `request.auth.uid == userId` in your rules

### "Firebase App named '[DEFAULT]' already exists"
- This usually means Firebase is initialized multiple times
- Check that `firebase.ts` is only imported once

## Next Steps

1. **Add Authentication UI**: Create login/signup components
2. **Wire up Library Page**: Update `Library.tsx` to use Firebase instead of backend API
3. **Add Save Button**: Update `Lab.tsx` to save molecules to Firebase
4. **Add Mini Viewer**: Use `MiniMoleculeViewer` component in molecule cards

## Resources

- [Firebase Documentation](https://firebase.google.com/docs)
- [Firestore Documentation](https://firebase.google.com/docs/firestore)
- [Firebase Auth Documentation](https://firebase.google.com/docs/auth)

