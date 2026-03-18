import type { ThreeEvent } from '@react-three/fiber'
import type { useLabStore } from '../store/labStore'
import type { ToolName } from '../types/molecule'

export type ToolPointerEvent = ThreeEvent<PointerEvent>
export type LabStoreState = ReturnType<typeof useLabStore.getState>

export interface Tool {
  name: ToolName
  onPointerDown?: (ev: ToolPointerEvent, store: LabStoreState) => void
  onPointerMove?: (ev: ToolPointerEvent, store: LabStoreState) => void
  onPointerUp?: (ev: ToolPointerEvent, store: LabStoreState) => void
}
