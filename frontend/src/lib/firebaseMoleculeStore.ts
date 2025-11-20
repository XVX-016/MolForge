import { db } from '../firebase';
import {
  collection,
  doc,
  setDoc,
  getDocs,
  deleteDoc,
  query,
  orderBy,
  where,
  Timestamp,
} from 'firebase/firestore';

export interface FirebaseMolecule {
  id?: string;
  name: string;
  smiles?: string;
  formula?: string;
  json_graph?: string;
  properties?: string;
  thumbnail_b64?: string;
  userId: string;
  createdAt: number;
}

/**
 * Save a molecule to Firestore for a specific user
 */
export async function saveMolecule(userId: string, molecule: Omit<FirebaseMolecule, 'id' | 'userId' | 'createdAt'>): Promise<string> {
  const ref = doc(collection(db, 'users', userId, 'molecules'));
  const data: Omit<FirebaseMolecule, 'id'> = {
    ...molecule,
    userId,
    createdAt: Date.now(),
  };
  await setDoc(ref, data);
  return ref.id;
}

/**
 * List all molecules for a specific user
 */
export async function listMolecules(userId: string): Promise<FirebaseMolecule[]> {
  const q = query(
    collection(db, 'users', userId, 'molecules'),
    orderBy('createdAt', 'desc')
  );
  const snap = await getDocs(q);
  return snap.docs.map((d) => ({
    id: d.id,
    ...d.data(),
  })) as FirebaseMolecule[];
}

/**
 * Delete a molecule for a specific user
 */
export async function deleteMolecule(userId: string, moleculeId: string): Promise<void> {
  await deleteDoc(doc(db, 'users', userId, 'molecules', moleculeId));
}

/**
 * Get a molecule by ID for a specific user
 */
export async function getMolecule(userId: string, moleculeId: string): Promise<FirebaseMolecule | null> {
  const docRef = doc(db, 'users', userId, 'molecules', moleculeId);
  const docSnap = await getDocs(collection(db, 'users', userId, 'molecules'));
  const found = docSnap.docs.find((d) => d.id === moleculeId);
  if (!found) return null;
  return {
    id: found.id,
    ...found.data(),
  } as FirebaseMolecule;
}

/**
 * Search molecules by name or SMILES for a specific user
 */
export async function searchMolecules(
  userId: string,
  searchQuery: string
): Promise<FirebaseMolecule[]> {
  const allMolecules = await listMolecules(userId);
  const query = searchQuery.toLowerCase().trim();
  if (!query) return allMolecules;
  
  return allMolecules.filter(
    (mol) =>
      mol.name.toLowerCase().includes(query) ||
      (mol.smiles || '').toLowerCase().includes(query) ||
      (mol.formula || '').toLowerCase().includes(query)
  );
}

