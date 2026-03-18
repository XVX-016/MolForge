import { beforeEach, describe, expect, it, vi } from 'vitest'
import { fireEvent, screen } from '@testing-library/react'
import type { TemplateData } from '../../../kernel/templateLoader'
import { renderWithProviders } from '../../../tests/test-utils'
import { TemplatePanel } from '../TemplatePanel'

vi.mock('../TemplatePanel.css', () => ({}))

const mockToggleTemplatePanel = vi.fn()
const mockStartTemplateDrag = vi.fn()
const mockStopTemplateDrag = vi.fn()

vi.mock('../../../store/uiStore', () => ({
  useUIStore: () => ({
    templatePanelExpanded: true,
    toggleTemplatePanel: mockToggleTemplatePanel,
  }),
}))

vi.mock('../../../store/templateTool.store', () => ({
  useTemplateToolStore: () => ({
    startTemplateDrag: mockStartTemplateDrag,
    stopTemplateDrag: mockStopTemplateDrag,
  }),
}))

const mockGetTemplates = vi.fn(() => [
  { id: 'water', name: 'Water (H2O)', data: { atoms: [], bonds: [] } },
  { id: 'benzene', name: 'Benzene (C6H6)', data: { atoms: [], bonds: [] } },
  { id: 'methane', name: 'Methane (CH4)', data: { atoms: [], bonds: [] } },
  { id: 'ethanol', name: 'Ethanol (C2H5OH)', data: { atoms: [], bonds: [] } },
])

const mockLoadTemplate = vi.fn((): TemplateData => ({ atoms: [], bonds: [] }))
const mockPlaceTemplate = vi.fn()

vi.mock('../../../kernel/templateLoader', () => ({
  getTemplates: () => mockGetTemplates(),
  loadTemplate: (templateId: string) => mockLoadTemplate(templateId),
  placeTemplate: (template: TemplateData, offset: { x: number; y: number; z: number }) => mockPlaceTemplate(template, offset),
}))

describe('TemplatePanel', () => {
  beforeEach(() => {
    vi.clearAllMocks()
  })

  it('renders template panel with header', () => {
    renderWithProviders(<TemplatePanel />)
    expect(screen.getByText('Templates')).toBeInTheDocument()
  })

  it('renders categories when expanded', () => {
    renderWithProviders(<TemplatePanel />)
    expect(screen.getByText('Rings')).toBeInTheDocument()
    expect(screen.getByText('Functional Groups')).toBeInTheDocument()
    expect(screen.getByText('Core Atoms')).toBeInTheDocument()
  })

  it('filters templates by search query', () => {
    renderWithProviders(<TemplatePanel />)
    const searchInput = screen.getByPlaceholderText('Search templates...')

    fireEvent.change(searchInput, { target: { value: 'water' } })

    expect(screen.getByText('Water (H2O)')).toBeInTheDocument()
  })

  it('calls placeTemplate when template is clicked', () => {
    renderWithProviders(<TemplatePanel />)

    const waterTemplate = screen.getByText('Water (H2O)')
    fireEvent.click(waterTemplate.closest('.template-item') ?? waterTemplate)

    expect(mockPlaceTemplate).toHaveBeenCalled()
  })

  it('handles drag start events', () => {
    renderWithProviders(<TemplatePanel />)

    const templateItem = screen.getByText('Water (H2O)').closest('[draggable]')
    if (templateItem) {
      fireEvent.dragStart(templateItem)
      expect(mockStartTemplateDrag).toHaveBeenCalled()
    }
  })

  it('handles drag end events', () => {
    renderWithProviders(<TemplatePanel />)

    const templateItem = screen.getByText('Water (H2O)').closest('[draggable]')
    if (templateItem) {
      fireEvent.dragEnd(templateItem)
      expect(mockStopTemplateDrag).toHaveBeenCalled()
    }
  })

  it('toggles panel expansion', () => {
    renderWithProviders(<TemplatePanel />)

    const toggleButton = screen.getByLabelText('Collapse')
    fireEvent.click(toggleButton)

    expect(mockToggleTemplatePanel).toHaveBeenCalled()
  })
})
