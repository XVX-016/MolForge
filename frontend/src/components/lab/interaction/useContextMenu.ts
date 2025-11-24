import { useState, useEffect, useRef } from 'react'
import { useLabStore } from '../../../store/labStore'

interface ContextMenuState {
  visible: boolean
  x: number
  y: number
  atomId: string | null
  bondId: string | null
}

/**
 * Hook for context menu handling
 */
export function useContextMenu() {
  const [menuState, setMenuState] = useState<ContextMenuState>({
    visible: false,
    x: 0,
    y: 0,
    atomId: null,
    bondId: null,
  })
  const store = useLabStore.getState()

  useEffect(() => {
    const handleContextMenu = (e: MouseEvent) => {
      e.preventDefault()
      
      // Check if right-clicked on atom or bond (would need to be passed from 3D scene)
      // For now, we'll handle it via the ToolHandler
    }

    const handleClick = () => {
      setMenuState(prev => ({ ...prev, visible: false }))
    }

    window.addEventListener('contextmenu', handleContextMenu)
    window.addEventListener('click', handleClick)

    return () => {
      window.removeEventListener('contextmenu', handleContextMenu)
      window.removeEventListener('click', handleClick)
    }
  }, [])

  const showMenu = (x: number, y: number, atomId?: string, bondId?: string) => {
    setMenuState({
      visible: true,
      x,
      y,
      atomId: atomId || null,
      bondId: bondId || null,
    })
  }

  const hideMenu = () => {
    setMenuState(prev => ({ ...prev, visible: false }))
  }

  const deleteAtom = (atomId: string) => {
    store.deleteAtom(atomId)
    hideMenu()
  }

  const deleteBond = (bondId: string) => {
    store.deleteBond(bondId)
    hideMenu()
  }

  return {
    menuState,
    showMenu,
    hideMenu,
    deleteAtom,
    deleteBond,
  }
}

