import React, { useState, useEffect } from 'react';

interface SearchBarProps {
  value?: string;
  onChange: (query: string) => void;
  placeholder?: string;
  debounceMs?: number;
}

/**
 * SearchBar - Debounced search input component
 * 
 * Automatically debounces input changes to avoid excessive API calls
 */
export default function SearchBar({
  value = '',
  onChange,
  placeholder = 'Search name, formula or SMILES...',
  debounceMs = 350,
}: SearchBarProps) {
  const [q, setQ] = useState(value);

  // Sync with external value changes
  useEffect(() => {
    setQ(value);
  }, [value]);

  // Debounce onChange callback
  useEffect(() => {
    const timer = setTimeout(() => {
      onChange(q.trim());
    }, debounceMs);

    return () => clearTimeout(timer);
  }, [q, onChange, debounceMs]);

  return (
    <div className="flex gap-3 items-center">
      <input
        value={q}
        onChange={(e) => setQ(e.target.value)}
        placeholder={placeholder}
        className="flex-1 rounded-full border border-gray-200 px-4 py-3 shadow-sm focus:outline-none focus:ring-2 focus:ring-darkGrey/20 focus:border-darkGrey placeholder:text-midGrey"
      />
      <button
        type="button"
        onClick={() => onChange(q.trim())}
        className="px-4 py-2 rounded-full bg-black text-white hover:bg-darkGrey transition-colors"
      >
        Search
      </button>
    </div>
  );
}

