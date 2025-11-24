import React from 'react'
import { useLabStore } from '../../store/labStore'
import AnimatedAtom from './AnimatedAtom'
import AnimatedBond from './AnimatedBond'
import BondMeshEnhanced from './viewer/BondMeshEnhanced'
import ToolHandler from './ToolHandler'
import { useDragging } from './interaction/useDragging'

export default function MoleculeScene(){
  const mol = useLabStore(s => s.molecule)
  const { startDrag } = useDragging()

  // Safety check
  if (!mol || !mol.atoms || !mol.bonds) {
    return null
  }

  return (
    <group>
      {/* render bonds first so atoms sit on top */}
      {mol.bonds.map((bond, i) => {
        if (!bond || !bond.atom1 || !bond.atom2) return null
        const a1 = mol.atoms.find(a => a && a.id === bond.atom1)
        const a2 = mol.atoms.find(a => a && a.id === bond.atom2)
        if (!a1 || !a2) return null
        return <BondMeshEnhanced key={bond.id || `bond-${i}`} bond={bond} atom1={a1} atom2={a2} index={i} />
      })}
      {mol.atoms.map((atom, i) => {
        if (!atom || !atom.id) return null
        return <AnimatedAtom key={atom.id} atom={atom} index={i} onDragStart={startDrag} />
      })}
      <ToolHandler />
    </group>
  )
}

