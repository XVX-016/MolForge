import React from 'react'
import { Atom, CloudUpload, Download, Link2, MousePointer2, Redo2, Trash2, Undo2, Upload } from 'lucide-react'
import { LibraryAPI } from '../../../api/library'
import { useLabStore } from '../../../store/labStore'
import type { Molecule, ToolName } from '../../../types/molecule'

type ActionId = ToolName | 'download' | 'upload' | 'save-library' | 'undo' | 'redo'

interface UploadFileTarget extends EventTarget {
  files?: FileList | null
}

function isMolecule(value: unknown): value is Molecule {
  if (!value || typeof value !== 'object') return false
  const candidate = value as Partial<Molecule>
  return Array.isArray(candidate.atoms) && Array.isArray(candidate.bonds)
}

export default function TopActionBar() {
  const { currentTool, setTool, molecule, loadMolecule } = useLabStore()

  const handleAction = async (id: ActionId) => {
    if (id === 'select' || id === 'add_atom' || id === 'bond' || id === 'delete') {
      setTool(id)
      return
    }

    if (id === 'download') {
      const dataStr = `data:text/json;charset=utf-8,${encodeURIComponent(JSON.stringify(molecule))}`
      const downloadAnchorNode = document.createElement('a')
      downloadAnchorNode.setAttribute('href', dataStr)
      downloadAnchorNode.setAttribute('download', `${molecule.metadata?.name || 'molecule'}.json`)
      document.body.appendChild(downloadAnchorNode)
      downloadAnchorNode.click()
      downloadAnchorNode.remove()
      return
    }

    if (id === 'upload') {
      const input = document.createElement('input')
      input.type = 'file'
      input.accept = '.json'
      input.onchange = (event: Event) => {
        const file = (event.target as UploadFileTarget).files?.[0]
        if (!file) return

        const reader = new FileReader()
        reader.onload = (readerEvent) => {
          try {
            const raw = readerEvent.target?.result
            if (typeof raw !== 'string') {
              throw new Error('Invalid file contents')
            }

            const parsed = JSON.parse(raw)
            if (!isMolecule(parsed)) {
              throw new Error('Invalid molecule JSON')
            }

            loadMolecule(parsed)
          } catch {
            alert('Invalid JSON file')
          }
        }
        reader.readAsText(file)
      }
      input.click()
      return
    }

    if (id === 'save-library') {
      const name = prompt('Enter molecule name', 'New Molecule')
      if (!name) return

      try {
        await LibraryAPI.upload({
          name,
          json_graph: { atoms: molecule.atoms, bonds: molecule.bonds },
        })
        alert('Saved to Library!')
      } catch (error) {
        console.error(error)
        alert('Failed to save.')
      }
      return
    }

    if (id === 'undo') {
      useLabStore.getState().undo()
      return
    }

    if (id === 'redo') {
      useLabStore.getState().redo()
    }
  }

  const actions: Array<{ id: ActionId; label: string; icon: React.ComponentType<{ size?: number; strokeWidth?: number; className?: string }> }> = [
    { id: 'select', label: 'Select', icon: MousePointer2 },
    { id: 'add_atom', label: 'Add Atom', icon: Atom },
    { id: 'bond', label: 'Add Bond', icon: Link2 },
    { id: 'delete', label: 'Delete', icon: Trash2 },
    { id: 'undo', label: 'Undo', icon: Undo2 },
    { id: 'redo', label: 'Redo', icon: Redo2 },
    { id: 'download', label: 'Download JSON', icon: Download },
    { id: 'upload', label: 'Upload JSON', icon: Upload },
    { id: 'save-library', label: 'Save to Cloud', icon: CloudUpload },
  ]

  return (
    <div className="w-full h-full flex items-center justify-center bg-white border-b border-gray-100 px-4">
      <div className="flex items-center gap-1 bg-gray-50/50 p-1.5 rounded-xl border border-gray-100">
        {actions.map((action) => {
          const Icon = action.icon
          const isTool = ['select', 'add_atom', 'bond', 'delete'].includes(action.id)
          const isActive = isTool && currentTool === action.id

          return (
            <button
              key={action.id}
              onClick={() => void handleAction(action.id)}
              className={`
                w-9 h-9 rounded-lg flex items-center justify-center transition-all duration-200
                ${isActive ? 'bg-white text-blue-600 shadow-sm ring-1 ring-gray-200' : 'text-gray-500 hover:bg-gray-100 hover:text-gray-900'}
              `}
              title={action.label}
            >
              <Icon size={18} strokeWidth={2} />
            </button>
          )
        })}
      </div>
    </div>
  )
}
