/**
 * Seed Library Page
 * Utility page to populate the library with sample molecules
 */
import React, { useState, useEffect } from 'react';
import { motion } from 'framer-motion';
import { supabase } from '../supabase';
import { seedMolecules, SAMPLE_MOLECULES } from '../utils/seedMolecules';

export default function SeedLibrary() {
  const [userId, setUserId] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);
  const [progress, setProgress] = useState({ current: 0, total: SAMPLE_MOLECULES.length });
  const [result, setResult] = useState<{ success: number; failed: number } | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (!supabase) return;

    // Check current session
    supabase.auth.getSession().then(({ data: { session } }) => {
      if (session?.user) {
        setUserId(session.user.id);
      } else {
        setUserId(null);
      }
    });

    // Listen for auth changes
    const {
      data: { subscription },
    } = supabase.auth.onAuthStateChange((_event, session) => {
      if (session?.user) {
        setUserId(session.user.id);
      } else {
        setUserId(null);
      }
    });

    return () => subscription.unsubscribe();
  }, []);

  const handleSeed = async () => {
    if (!userId) {
      setError('Please sign in to seed molecules');
      return;
    }

    if (!confirm(`This will add ${SAMPLE_MOLECULES.length} sample molecules to your library. Continue?`)) {
      return;
    }

    setLoading(true);
    setError(null);
    setResult(null);
    setProgress({ current: 0, total: SAMPLE_MOLECULES.length });

    try {
      const result = await seedMolecules(userId, (current, total) => {
        setProgress({ current, total });
      });
      setResult(result);
    } catch (err: any) {
      setError(err.message || 'Failed to seed molecules');
      console.error('Seeding error:', err);
    } finally {
      setLoading(false);
    }
  };

  if (!userId) {
    return (
      <motion.div
        initial={{ opacity: 0 }}
        animate={{ opacity: 1 }}
        className="p-8 space-y-6"
      >
        <div className="text-center py-12">
          <h2 className="text-2xl font-bold text-black mb-2">Authentication Required</h2>
          <p className="text-midGrey">Please sign in to seed your molecule library.</p>
        </div>
      </motion.div>
    );
  }

  return (
    <motion.div
      initial={{ opacity: 0 }}
      animate={{ opacity: 1 }}
      className="p-8 space-y-6 max-w-2xl mx-auto"
    >
      <header>
        <h1 className="text-3xl font-bold text-black mb-2">Seed Library</h1>
        <p className="text-darkGrey">
          Populate your library with {SAMPLE_MOLECULES.length} sample molecules including:
        </p>
        <ul className="list-disc list-inside text-midGrey mt-2 space-y-1">
          <li>Aspirin, Caffeine, Glucose, Ethanol</li>
          <li>Paracetamol, Ibuprofen, Vitamin C</li>
          <li>Nicotine, Serotonin, Dopamine</li>
        </ul>
      </header>

      <div className="bg-white rounded-xl shadow-neon border border-lightGrey p-6 space-y-4">
        <div>
          <h2 className="text-lg font-semibold text-black mb-2">Sample Molecules</h2>
          <div className="grid grid-cols-2 md:grid-cols-3 gap-2 text-sm">
            {SAMPLE_MOLECULES.map((mol) => (
              <div key={mol.name} className="text-darkGrey">
                {mol.name} ({mol.formula})
              </div>
            ))}
          </div>
        </div>

        {loading && (
          <div className="space-y-2">
            <div className="flex items-center justify-between text-sm">
              <span className="text-darkGrey">Generating molecules...</span>
              <span className="text-midGrey">
                {progress.current} / {progress.total}
              </span>
            </div>
            <div className="w-full bg-gray-200 rounded-full h-2">
              <motion.div
                className="bg-black h-2 rounded-full"
                initial={{ width: 0 }}
                animate={{ width: `${(progress.current / progress.total) * 100}%` }}
                transition={{ duration: 0.3 }}
              />
            </div>
          </div>
        )}

        {error && (
          <div className="bg-red-50 border border-red-200 rounded-lg p-4 text-red-800 text-sm">
            {error}
          </div>
        )}

        {result && (
          <div className="bg-green-50 border border-green-200 rounded-lg p-4 space-y-2">
            <div className="text-green-800 font-semibold">Seeding Complete!</div>
            <div className="text-green-700 text-sm">
              <div>✓ Successfully added: {result.success} molecules</div>
              {result.failed > 0 && <div>✗ Failed: {result.failed} molecules</div>}
            </div>
            <div className="text-green-600 text-xs mt-2">
              Visit the Library page to see your new molecules.
            </div>
          </div>
        )}

        <button
          onClick={handleSeed}
          disabled={loading}
          className="w-full px-4 py-3 rounded-lg bg-black text-white font-medium hover:bg-darkGrey disabled:opacity-50 disabled:cursor-not-allowed transition-all"
        >
          {loading ? 'Generating Molecules...' : `Seed Library (${SAMPLE_MOLECULES.length} molecules)`}
        </button>
      </div>

      <div className="bg-blue-50 border border-blue-200 rounded-lg p-4 text-sm text-blue-800">
        <strong>Note:</strong> This will generate 3D molfiles and thumbnails for each molecule.
        This may take a few minutes depending on your backend API speed.
      </div>
    </motion.div>
  );
}

