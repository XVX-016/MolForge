import { useEffect } from 'react'
import { useLabStore } from '../../../store/labStore'

/**
 * Hook for keyboard shortcuts (undo/redo)
 */
export function useKeyboardShortcuts() {
  const undo = useLabStore(s => s.undo)
  const redo = useLabStore(s => s.redo)

  useEffect(() => {
    const handleKeyDown = (event: KeyboardEvent) => {
      // Ctrl+Z or Cmd+Z for undo
      if ((event.ctrlKey || event.metaKey) && event.key === 'z' && !event.shiftKey) {
        event.preventDefault()
        undo()
      }
      // Ctrl+Shift+Z or Cmd+Shift+Z for redo
      if ((event.ctrlKey || event.metaKey) && event.key === 'z' && event.shiftKey) {
        event.preventDefault()
        redo()
      }
    }

    window.addEventListener('keydown', handleKeyDown)
    return () => {
      window.removeEventListener('keydown', handleKeyDown)
    }
  }, [undo, redo])
}

