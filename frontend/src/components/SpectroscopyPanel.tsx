import { useMemo, useState } from 'react'
import { useMoleculeStore } from '../store/moleculeStore'
import IRGraph from './spectroscopy/IRGraph'
import NMRGraph from './spectroscopy/NMRGraph'
import MassGraph from './spectroscopy/MassGraph'

type TabId = 'ir' | 'nmr' | 'mass'

const tabs: { id: TabId; label: string }[] = [
  { id: 'ir', label: 'IR Spectrum' },
  { id: 'nmr', label: 'NMR' },
  { id: 'mass', label: 'Mass Spec' },
]

export default function SpectroscopyPanel() {
  const [activeTab, setActiveTab] = useState<TabId>('ir')
  const spectroscopy = useMoleculeStore((s) => s.spectroscopy)
  const runSpectroscopy = useMoleculeStore((s) => s.runSpectroscopy)
  const setHighlightedAtoms = useMoleculeStore((s) => s.setHighlightedAtoms)
  const molecule = useMoleculeStore((s) => s.currentMolecule)

  const ready = Boolean(molecule && molecule.atoms.size > 0)
  const irPeaks = useMemo(() => spectroscopy.ir?.peaks ?? [], [spectroscopy.ir])
  const protonPeaks = useMemo(() => spectroscopy.nmr?.protons ?? [], [spectroscopy.nmr])
  const carbonPeaks = useMemo(() => spectroscopy.nmr?.carbon13 ?? [], [spectroscopy.nmr])
  const massPeaks = useMemo(() => spectroscopy.mass?.peaks ?? [], [spectroscopy.mass])

  const functionalGroups = useMemo(() => {
    const groups = new Map<string, number>()
    irPeaks.forEach((peak) => {
      groups.set(peak.bondType, peak.frequency)
    })
    return Array.from(groups.entries()).slice(0, 4)
  }, [irPeaks])

  const statusLabel = (() => {
    if (!ready) return 'Add atoms to calculate spectra'
    if (spectroscopy.status === 'calculating') return 'Calculating...'
    if (spectroscopy.status === 'complete') return 'Spectra ready'
    if (spectroscopy.status === 'error') return spectroscopy.error ?? 'Unable to calculate spectra'
    return 'Idle'
  })()

  const handlePeakHover = (atomIds?: string[]) => {
    setHighlightedAtoms(atomIds ?? [])
  }

  const formattedTimestamp = spectroscopy.timestamp
    ? new Intl.DateTimeFormat('en', {
      hour: 'numeric',
      minute: 'numeric',
      second: 'numeric',
    }).format(new Date(spectroscopy.timestamp))
    : null

  return (
    <div className="space-y-4">
      <div className="flex items-center justify-between gap-2">
        <div>
          <p className="text-xs uppercase tracking-[0.3em] text-midGrey">Spectroscopy</p>
          <p className="text-sm text-darkGrey">{statusLabel}</p>
          {formattedTimestamp && (
            <p className="text-[11px] text-midGrey mt-0.5">Updated {formattedTimestamp}</p>
          )}
        </div>
        <button
          onClick={runSpectroscopy}
          className="btn-secondary text-xs px-3 py-1.5"
          disabled={!ready || spectroscopy.status === 'calculating'}
        >
          {spectroscopy.status === 'calculating' ? 'Working…' : 'Recalculate'}
        </button>
      </div>

      {!ready ? (
        <p className="text-sm text-midGrey">Build a molecule to preview spectra.</p>
      ) : (
        <>
          <div className="flex gap-2 border-b border-lightGrey">
            {tabs.map((tab) => (
              <button
                key={tab.id}
                onClick={() => setActiveTab(tab.id)}
                className={`px-3 py-2 text-xs font-semibold transition ${activeTab === tab.id
                    ? 'text-black border-b-2 border-black'
                    : 'text-midGrey hover:text-black'
                  }`}
              >
                {tab.label}
              </button>
            ))}
          </div>

          <div className="space-y-4">
            {activeTab === 'ir' && irPeaks.length > 0 && (
              <>
                <IRGraph
                  peaks={irPeaks.map((peak) => ({
                    wavenumber: peak.frequency,
                    intensity: peak.intensity > 0.75 ? 'Strong' : peak.intensity > 0.5 ? 'Medium' : 'Weak',
                    group: peak.bondType,
                    atoms: peak.atoms,
                  }))}
                  onPeakHover={(peak) => handlePeakHover(peak?.atoms)}
                />
                <PeakTable
                  headers={['Frequency (cm⁻¹)', 'Bond', 'Annotation']}
                  rows={irPeaks.map((peak) => [
                    peak.frequency.toFixed(0),
                    peak.bondType,
                    peak.annotation ?? '—',
                  ])}
                />
              </>
            )}

            {activeTab === 'nmr' && (protonPeaks.length > 0 || carbonPeaks.length > 0) && (
              <div className="space-y-4">
                {protonPeaks.length > 0 && (
                  <div className="space-y-2">
                    <p className="text-xs uppercase tracking-[0.2em] text-midGrey">¹H NMR</p>
                    <NMRGraph
                      peaks={protonPeaks.map((peak) => ({
                        shift: parseFloat(peak.chemicalShift.toFixed(2)),
                        multiplicity: peak.multiplicity,
                        integration: peak.integration,
                        atoms: peak.atoms,
                      }))}
                      type="1H"
                      onPeakHover={(peak) => handlePeakHover(peak?.atoms)}
                    />
                  </div>
                )}
                {carbonPeaks.length > 0 && (
                  <div className="space-y-2">
                    <p className="text-xs uppercase tracking-[0.2em] text-midGrey">¹³C NMR</p>
                    <NMRGraph
                      peaks={carbonPeaks.map((peak) => ({
                        shift: parseFloat(peak.chemicalShift.toFixed(1)),
                        multiplicity: 's',
                        integration: 1,
                        atoms: peak.atoms,
                      }))}
                      type="13C"
                      onPeakHover={(peak) => handlePeakHover(peak?.atoms)}
                    />
                  </div>
                )}
              </div>
            )}

            {activeTab === 'mass' && massPeaks.length > 0 && spectroscopy.mass && (
              <>
                <MassGraph
                  peaks={massPeaks.map((peak) => ({
                    'm/z': peak.mz,
                    intensity: peak.intensity,
                    fragment_type: peak.fragment,
                    fragment_atoms: peak.atoms,
                  }))}
                  molecularIon={spectroscopy.mass.parentMass}
                  onPeakHover={(peak) => handlePeakHover(peak?.fragment_atoms)}
                />
                <PeakTable
                  headers={['m/z', 'Intensity', 'Fragment']}
                  rows={massPeaks.map((peak) => [
                    peak.mz.toFixed(1),
                    `${Math.round(peak.intensity * 100)}%`,
                    peak.fragment,
                  ])}
                />
              </>
            )}
          </div>

          {functionalGroups.length > 0 && (
            <div className="pt-2 border-t border-lightGrey">
              <p className="text-xs uppercase tracking-[0.2em] text-midGrey mb-2">
                Functional highlights
              </p>
              <div className="flex flex-wrap gap-2">
                {functionalGroups.map(([bond, freq]) => (
                  <span
                    key={bond}
                    className="text-xs px-2 py-1 rounded-full bg-offwhite border border-lightGrey"
                  >
                    {bond} · {freq.toFixed(0)} cm⁻¹
                  </span>
                ))}
              </div>
            </div>
          )}
        </>
      )}
    </div>
  )
}

interface PeakTableProps {
  headers: string[]
  rows: (string | number)[][]
}

function PeakTable({ headers, rows }: PeakTableProps) {
  if (rows.length === 0) return null
  return (
    <div className="overflow-x-auto border border-lightGrey rounded-lg">
      <table className="min-w-full text-sm">
        <thead>
          <tr className="bg-offwhite text-left text-xs uppercase tracking-wide text-midGrey">
            {headers.map((header) => (
              <th key={header} className="px-3 py-2 font-semibold">
                {header}
              </th>
            ))}
          </tr>
        </thead>
        <tbody>
          {rows.slice(0, 6).map((row, index) => (
            <tr
              key={index}
              className="border-t border-lightGrey text-sm text-darkGrey odd:bg-white even:bg-offwhite/60"
            >
              {row.map((cell, cellIndex) => (
                <td key={cellIndex} className="px-3 py-2">
                  {cell}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  )
}
