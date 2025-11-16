/// <reference types="vitest/config" />

import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';
import { configDefaults } from 'vitest/config';
import path from 'path';
import { fileURLToPath } from 'url';

const __dirname = path.dirname(fileURLToPath(import.meta.url));

export default defineConfig({
  // @ts-expect-error - Vitest uses its own Vite instance, but plugins are compatible at runtime
  plugins: [
    react(),
  ],

  resolve: {
    alias: {
      '@biosynth/engine': path.resolve(__dirname, '../packages/engine/src/index.ts'),
    },
  },

  test: {
    globals: true,
    environment: 'jsdom',
    setupFiles: ['./src/tests/vitest.setup.ts'],
    exclude: [...configDefaults.exclude],
  },
});
