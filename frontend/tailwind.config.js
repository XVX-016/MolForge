/** @type {import('tailwindcss').Config} */
export default {
  content: [
    "./index.html",
    "./src/**/*.{js,ts,jsx,tsx}",
  ],
  theme: {
    extend: {
      colors: {
        // Primary Colors - MUST use CSS variables
        chrome: 'var(--color-chrome)',
        ivory: 'var(--color-ivory)',
        spaceGrey: 'var(--color-spaceGrey)',
        ionBlack: 'var(--color-ionBlack)',
        frostedGlass: 'var(--color-frostedGlass)',
        
        // Accents
        neonCyan: 'var(--color-neonCyan)',
        plasmaTeal: 'var(--color-plasmaTeal)',
        violetEdge: 'var(--color-violetEdge)',
        
        // Legacy aliases (deprecated - use direct colors)
        aluminum: {
          light: 'var(--color-ivory)',
          DEFAULT: 'var(--color-chrome)',
          dark: 'var(--color-spaceGrey)',
        },
        panel: {
          DEFAULT: 'var(--color-frostedGlass)',
          soft: 'rgba(255,255,255,0.03)',
        },
        accent: {
          cyan: 'var(--color-neonCyan)',
          teal: 'var(--color-plasmaTeal)',
          violet: 'var(--color-violetEdge)',
        },
        text: {
          primary: 'var(--color-ivory)',
          secondary: 'var(--color-chrome)',
          tertiary: 'var(--color-spaceGrey)',
        },
      },
      backgroundImage: {
        'chrome-gradient': 'var(--gradient-chrome)',
        'ivory-gradient': 'var(--gradient-ivory)',
        'plasma-neon': 'var(--gradient-plasmaNeon)',
      },
      borderRadius: {
        xl: '22px',
        lg: '18px',
      },
      boxShadow: {
        chrome: 'var(--shadow-chrome)',
        glass: 'var(--shadow-glass)',
        neon: 'var(--shadow-neon)',
        'neon-sm': 'var(--shadow-neon-hover)',
      },
    },
  },
  plugins: [],
}
