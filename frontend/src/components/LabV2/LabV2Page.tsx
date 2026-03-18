import { useEffect } from 'react'
import { useLocation, useSearchParams } from 'react-router-dom'

import BottomTemplateBar from './layout/BottomTemplateBar'
import LabCanvas from './LabCanvas'
import LeftToolDock from './layout/LeftToolDock'
import RightInspector from './layout/RightInspector'

import { getPublicMolecule, type PublicMolecule } from '../../lib/publicMoleculeStore'
import { getUserMolecule, type UserMolecule } from '../../lib/userMoleculeStore'
import { supabase } from '../../supabase'
import { useLabStore } from '../../store/labStore'
import { parseMoleculeFromSupabase } from '../../utils/moleculeParser'

type NavigationState = {
  source?: 'user' | 'public'
  moleculeId?: string | number
  molfile?: string | null
  name?: string
  smiles?: string | null
  formula?: string | null
}

type LoadableMoleculeData = {
  id?: string | number
  name?: string
  smiles?: string | null
  formula?: string | null
  molfile?: string | null
  json_graph?: string | null
  source?: 'user' | 'public'
}

function withJsonGraph<T extends { id?: string; name?: string; smiles?: string; formula?: string; molfile?: string }>(
  molecule: T,
  source: 'user' | 'public'
): LoadableMoleculeData {
  const jsonGraph = 'json_graph' in molecule && typeof molecule.json_graph === 'string' ? molecule.json_graph : null

  return {
    id: molecule.id,
    name: molecule.name,
    smiles: molecule.smiles,
    formula: molecule.formula,
    molfile: molecule.molfile,
    json_graph: jsonGraph,
    source,
  }
}

export default function LabV2Page() {
  const location = useLocation()
  const [searchParams] = useSearchParams()
  const loadMolecule = useLabStore((state) => state.loadMolecule)

  useEffect(() => {
    const loadFromNavigation = async () => {
      const state = location.state as NavigationState | undefined

      let source = state?.source
      let moleculeId = state?.moleculeId

      if (!source) {
        const qpSource = searchParams.get('source')
        if (qpSource === 'user' || qpSource === 'public') {
          source = qpSource
        }
      }

      if (!moleculeId) {
        const qpId = searchParams.get('id')
        if (qpId) {
          moleculeId = qpId
        }
      }

      if (!source || !moleculeId) {
        console.log('[LabV2Page] Missing source or moleculeId:', { source, moleculeId })
        return
      }

      console.log('[LabV2Page] Fetching molecule:', { source, moleculeId })

      try {
        let moleculeData: LoadableMoleculeData | null = null

        if (source === 'user') {
          if (!supabase) return

          const {
            data: { session },
          } = await supabase.auth.getSession()
          const userId = session?.user?.id
          if (!userId) return

          const userMol: UserMolecule | null = await getUserMolecule(userId, String(moleculeId))
          if (userMol) {
            moleculeData = withJsonGraph(userMol, 'user')
            console.log('[LabV2Page] Fetched user molecule data:', {
              id: moleculeData.id,
              name: moleculeData.name,
              hasMolfile: !!moleculeData.molfile,
              hasJsonGraph: !!moleculeData.json_graph,
            })
          } else {
            console.warn('[LabV2Page] User molecule not found:', { userId, moleculeId })
          }
        } else if (source === 'public') {
          const publicMol: PublicMolecule | null = await getPublicMolecule(String(moleculeId))
          if (publicMol) {
            moleculeData = withJsonGraph(publicMol, 'public')
            console.log('[LabV2Page] Fetched public molecule data:', {
              id: moleculeData.id,
              name: moleculeData.name,
              hasMolfile: !!moleculeData.molfile,
              hasJsonGraph: !!moleculeData.json_graph,
            })
          } else {
            console.warn('[LabV2Page] Public molecule not found:', { moleculeId })
          }
        }

        if (moleculeData) {
          if (!moleculeData.molfile && !moleculeData.json_graph && moleculeData.smiles?.trim()) {
            console.log('[LabV2Page] Generating molfile from SMILES:', moleculeData.smiles)
            try {
              const { convertSMILESToMolfile } = await import('../../lib/api')
              const result = await convertSMILESToMolfile(moleculeData.smiles.trim())
              if (result.molfile) {
                moleculeData.molfile = result.molfile
                console.log('[LabV2Page] Successfully generated molfile from SMILES via backend API')
              }
            } catch (error) {
              console.warn('[LabV2Page] Backend SMILES conversion failed, will use frontend parser fallback:', error)
            }
          }

          console.log('[LabV2Page] Parsing molecule from Supabase data')
          const molecule = parseMoleculeFromSupabase(moleculeData)
          console.log('[LabV2Page] Parsed molecule:', {
            id: molecule.id,
            atomCount: molecule.atoms.length,
            bondCount: molecule.bonds.length,
            name: molecule.metadata?.name,
          })
          console.log('[LabV2Page] Loading molecule into store')
          loadMolecule(molecule)
        } else if (state?.molfile) {
          console.log('[LabV2Page] Using fallback molfile from state')
          const molecule = parseMoleculeFromSupabase({
            id: moleculeId,
            name: state.name,
            smiles: state.smiles,
            formula: state.formula,
            molfile: state.molfile,
            source,
          })
          console.log('[LabV2Page] Parsed fallback molecule:', {
            id: molecule.id,
            atomCount: molecule.atoms.length,
            bondCount: molecule.bonds.length,
          })
          loadMolecule(molecule)
        } else {
          console.warn('[LabV2Page] No molecule data available to load')
        }
      } catch (error) {
        console.error('[LabV2Page] Failed to load molecule for Lab', error)
      }
    }

    void loadFromNavigation()
  }, [loadMolecule, location.state, searchParams])

  return (
    <div className="h-full w-full grid grid-rows-[1fr_72px] overflow-hidden bg-white">
      <div className="grid grid-cols-[72px_1fr] min-h-0 relative z-0">
        <div className="h-full z-10">
          <LeftToolDock />
        </div>

        <div className="relative w-full h-full bg-[#f8f9fa] overflow-hidden">
          <div className="absolute inset-0 z-0">
            <LabCanvas />
          </div>
          <RightInspector />
        </div>
      </div>

      <div className="z-20 border-t border-gray-100 bg-white">
        <BottomTemplateBar />
      </div>
    </div>
  )
}
