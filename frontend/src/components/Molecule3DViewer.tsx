import React, { useEffect, useRef } from 'react';

interface Molecule3DViewerProps {
  molfile?: string | null;
  smiles?: string | null;
  height?: number;
  backgroundColor?: string;
  autorotate?: boolean;
  className?: string;
  onLoaded?: () => void;
}

/**
 * Molecule3DViewer - Renders 3D molecules using 3Dmol.js with ball-and-stick representation
 * 
 * Supports both molfile (preferred) and SMILES input.
 * Falls back to thumbnail if 3D rendering fails.
 */
export default function Molecule3DViewer({
  molfile,
  smiles,
  height = 200,
  backgroundColor = 'white',
  autorotate = false,
  className = '',
  onLoaded,
}: Molecule3DViewerProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<any>(null);

  useEffect(() => {
    // Wait for 3Dmol to be loaded on window
    if (!window.$3Dmol) {
      console.warn('3Dmol.js not found. Make sure you included the CDN script in index.html.');
      return;
    }

    if (!containerRef.current) return;

    // Cleanup previous viewer if exists
    if (viewerRef.current) {
      try {
        viewerRef.current.clear();
      } catch (e) {
        // Ignore cleanup errors
      }
      viewerRef.current = null;
      if (containerRef.current) {
        containerRef.current.innerHTML = '';
      }
    }

    // Create new viewer
    try {
      const config = { backgroundColor };
      const viewer = window.$3Dmol.createViewer(containerRef.current, config);
      viewerRef.current = viewer;

      // Try to load molfile first (preferred)
      if (molfile) {
        try {
          viewer.addModel(molfile, 'mol');
          // Use ball-and-stick representation (stick with small radius + atom spheres)
          viewer.setStyle({}, {
            stick: { radius: 0.12, colorscheme: 'element' },
            sphere: { radius: 0.25, colorscheme: 'element' }
          });
          viewer.zoomTo();
          viewer.render();
          if (autorotate) {
            viewer.animate({ loop: true });
          }
          // Notify parent that viewer is loaded
          if (onLoaded) {
            // Small delay to ensure render is complete
            setTimeout(() => onLoaded(), 100);
          }
        } catch (err) {
          console.error('3Dmol render error with molfile:', err);
          // Fallback: try just sticks
          try {
            viewer.setStyle({}, { stick: { colorscheme: 'element' } });
            viewer.render();
            if (onLoaded) setTimeout(() => onLoaded(), 100);
          } catch (e) {
            console.error('3Dmol fallback render error:', e);
            if (onLoaded) onLoaded(); // Still notify even on error
          }
        }
      } else if (smiles) {
        // Try to load from SMILES (less reliable, but works for simple molecules)
        try {
          viewer.addModel(smiles, 'smi');
          viewer.setStyle({}, {
            stick: { radius: 0.12, colorscheme: 'element' },
            sphere: { radius: 0.25, colorscheme: 'element' }
          });
          viewer.zoomTo();
          viewer.render();
          if (autorotate) {
            viewer.animate({ loop: true });
          }
          if (onLoaded) {
            setTimeout(() => onLoaded(), 100);
          }
        } catch (err) {
          console.error('3Dmol render error with SMILES:', err);
          if (onLoaded) onLoaded(); // Still notify even on error
        }
      }

      // Handle window resize
      const handleResize = () => {
        if (viewerRef.current) {
          try {
            viewerRef.current.resize();
          } catch (e) {
            // Ignore resize errors
          }
        }
      };
      window.addEventListener('resize', handleResize);

      return () => {
        window.removeEventListener('resize', handleResize);
        if (viewerRef.current) {
          try {
            viewerRef.current.clear();
          } catch (e) {
            // Ignore cleanup errors
          }
          viewerRef.current = null;
        }
      };
    } catch (err) {
      console.error('Failed to create 3Dmol viewer:', err);
      if (onLoaded) onLoaded(); // Notify even on error
    }
  }, [molfile, smiles, backgroundColor, autorotate, onLoaded]);

  if (!molfile && !smiles) {
    return (
      <div
        className={`flex items-center justify-center bg-gray-50 rounded-md ${className}`}
        style={{ height }}
      >
        <div className="text-gray-400 text-sm">No 3D data</div>
      </div>
    );
  }

  return (
    <div
      ref={containerRef}
      className={`rounded-md overflow-hidden ${className}`}
      style={{ width: '100%', height }}
    />
  );
}

// Extend Window interface for 3Dmol
declare global {
  interface Window {
    $3Dmol?: {
      createViewer: (element: HTMLElement, config: { backgroundColor?: string }) => any;
    };
  }
}


