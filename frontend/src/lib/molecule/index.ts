/**
 * Molecule Editor Library
 * 
 * Centralized molecule editing, validation, and operations.
 * 
 * This is the single source of truth for molecule manipulation.
 */

// Core types
export type {
  Atom,
  Bond,
  MoleculeState,
  EditorTool,
  ValidationError,
  ValidationErrorCode,
  ValidationResult,
  ElementInfo,
  ExportOptions,
} from './types'

// Core classes
export { AtomImpl } from './Atom'
export { BondImpl } from './Bond'
export { Molecule } from './Molecule'

// Constants
export { ELEMENT_DATA, COMMON_ELEMENTS, DEFAULT_BOND_ORDER, VALID_BOND_ORDERS } from './constants'

// Re-export validation (will be created in Phase 5)
// export { validateMolecule } from './validation/Validator'

