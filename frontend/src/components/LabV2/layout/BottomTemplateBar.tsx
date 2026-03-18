import { Circle, Hexagon, Share2, Triangle } from 'lucide-react'
import { nanoid } from 'nanoid'
import { useLabStore } from '../../../store/labStore'
import type { Atom, Bond, MoleculeGraph } from '../../../types/molecule'
import { BENZENE, CHAIN_4, CYCLOHEXANE, CYCLOPROPANE } from '../../../utils/defaultMolecules'

const TEMPLATES = [
  { id: 'benzene', name: 'Benzene', data: BENZENE, icon: Hexagon },
  { id: 'cyclohexane', name: 'Cyclohexane', data: CYCLOHEXANE, icon: Circle },
  { id: 'chain', name: 'Chain', data: CHAIN_4, icon: Share2 },
  { id: 'triangle', name: 'Cyclopropane', data: CYCLOPROPANE, icon: Triangle },
]

interface LegacyAtom {
  id: string
  element: string
  x?: number
  y?: number
  z?: number
  position?: [number, number, number]
}

interface LegacyBond {
  from?: string
  to?: string
  a?: string
  b?: string
  order?: 1 | 1.5 | 2 | 3
}

type TemplateData = MoleculeGraph | { atoms: LegacyAtom[]; bonds: LegacyBond[] }

function normalizeTemplateData(data: TemplateData): MoleculeGraph {
  const atomMap = new Map<string, string>()
  const rawAtoms = data.atoms as Array<Atom | LegacyAtom>
  const rawBonds = data.bonds as Array<Bond | LegacyBond>

  const atoms: Atom[] = rawAtoms.map((atom) => {
    const newId = `atom-${nanoid()}`
    atomMap.set(atom.id, newId)

    const position = Array.isArray(atom.position)
      ? { x: atom.position[0], y: atom.position[1], z: atom.position[2] }
      : { x: atom.x ?? 0, y: atom.y ?? 0, z: atom.z ?? 0 }

    return {
      id: newId,
      element: atom.element,
      position,
      charge: 'charge' in atom && typeof atom.charge === 'number' ? atom.charge : 0,
    }
  })

  const bonds: Bond[] = rawBonds
    .map((bond) => {
      const from = atomMap.get('from' in bond && bond.from ? bond.from : bond.a ?? '')
      const to = atomMap.get('to' in bond && bond.to ? bond.to : bond.b ?? '')

      if (!from || !to) return null

      return {
        id: `bond-${nanoid()}`,
        from,
        to,
        order: bond.order ?? 1,
      }
    })
    .filter((bond): bond is Bond => bond !== null)

  return {
    id: `mol-${Date.now()}`,
    atoms,
    bonds,
    metadata: {},
  }
}

export default function BottomTemplateBar() {
  const { loadMolecule } = useLabStore()

  return (
    <div className="w-full h-full bg-white border-t border-gray-100 px-6 flex items-center gap-6 overflow-x-auto">
      <span className="text-[10px] uppercase font-bold text-gray-400 tracking-wider shrink-0">Templates</span>

      <div className="h-8 w-px bg-gray-200 shrink-0" />

      <div className="flex items-center gap-4">
        {TEMPLATES.map((template) => {
          const Icon = template.icon
          return (
            <button
              key={template.id}
              onClick={() => {
                loadMolecule(normalizeTemplateData(template.data))
              }}
              className="flex items-center gap-2 px-3 py-1.5 rounded-lg hover:bg-gray-50 text-gray-600 transition-colors border border-transparent hover:border-gray-200 group"
            >
              <Icon size={16} className="text-gray-400 group-hover:text-blue-500 transition-colors" />
              <span className="text-sm font-medium">{template.name}</span>
            </button>
          )
        })}
      </div>
    </div>
  )
}
