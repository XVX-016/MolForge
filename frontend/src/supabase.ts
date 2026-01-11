import { createClient, type SupabaseClient } from '@supabase/supabase-js';

// Supabase configuration from environment variables
const supabaseUrl = import.meta.env.VITE_SUPABASE_URL || import.meta.env.NEXT_PUBLIC_SUPABASE_URL;
const supabaseAnonKey = import.meta.env.VITE_SUPABASE_ANON_KEY || import.meta.env.NEXT_PUBLIC_SUPABASE_ANON_KEY;

// Debug: Log what Vite is actually reading (only in dev mode)
if (import.meta.env.DEV) {
  console.log('🔍 Debug - Environment variables from Vite:');
  console.log('  VITE_SUPABASE_URL:', supabaseUrl);
  console.log('  VITE_SUPABASE_ANON_KEY:', supabaseAnonKey ? 'Present (Hidden)' : 'Missing');

  if (!supabaseUrl || !supabaseAnonKey) {
    console.error("CRITICAL: Environment variables missing. Check .env file location.");
  }
}

// Validate Supabase configuration
function validateSupabaseConfig() {
  if (!supabaseUrl || !supabaseAnonKey) {
    console.warn(
      '⚠️ Supabase configuration incomplete. Missing:',
      !supabaseUrl ? 'VITE_SUPABASE_URL / NEXT_PUBLIC_SUPABASE_URL' : '',
      !supabaseAnonKey ? 'VITE_SUPABASE_ANON_KEY / NEXT_PUBLIC_SUPABASE_ANON_KEY' : '',
      '\nPlease create a .env file with your Supabase credentials.'
    );
    return false;
  }

  // Check for placeholder values
  const hasPlaceholders = supabaseUrl.includes('your-') ||
    supabaseAnonKey.includes('your-') ||
    supabaseUrl.includes('placeholder') ||
    supabaseAnonKey.includes('placeholder');

  if (hasPlaceholders) {
    console.warn('⚠️ Supabase configuration contains placeholder values.');
    console.warn('📝 Please update your .env file with actual Supabase credentials:');
    console.warn('   1. Go to https://supabase.com/dashboard');
    console.warn('   2. Select your project → Settings → API');
    console.warn('   3. Copy Project URL and anon/public key');
    console.warn('   4. Update frontend/.env with real values');
    console.warn('   5. Restart dev server');
    return false;
  }

  // Validate URL format
  if (!supabaseUrl.startsWith('https://')) {
    console.warn('⚠️ VITE_SUPABASE_URL appears invalid. Should start with https://');
    // return false; // Let's allow it to try anyway
  }

  /*
  // Strict domain check removed to allow self-hosted or other domains
  if (!supabaseUrl.includes('.supabase.co')) {
     console.warn('⚠️ VITE_SUPABASE_URL appears invalid. Should be: https://xxxxx.supabase.co');
  }
  */

  // Validate key format (anon keys start with eyJ, NOT sb_secret_)
  if (supabaseAnonKey.startsWith('sb_secret_')) {
    console.error('❌ ERROR: You are using a SECRET KEY instead of ANON KEY!');
    return false;
  }

  if (!supabaseAnonKey.startsWith('eyJ')) {
    console.warn('⚠️ VITE_SUPABASE_ANON_KEY appears invalid. Should start with "eyJ"');
    // return false; // Allow it anyway
  }

  return true;
}

let supabase: SupabaseClient | null = null;

try {
  const isValid = validateSupabaseConfig();

  if (isValid && supabaseUrl && supabaseAnonKey) {
    supabase = createClient(supabaseUrl, supabaseAnonKey, {
      auth: {
        persistSession: true,
        autoRefreshToken: true,
        detectSessionInUrl: true,
      }
    });
    console.log('✅ Supabase initialized successfully');
  } else {
    console.warn('⚠️ Supabase not initialized - configuration missing');
  }
} catch (error) {
  console.error('❌ Failed to initialize Supabase:', error);
}

// Export with null check
export { supabase };
export const isSupabaseConfigured = () => supabase !== null;

