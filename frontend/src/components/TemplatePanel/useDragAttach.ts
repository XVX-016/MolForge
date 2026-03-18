import { useCallback } from 'react'
import type { TemplateData } from '../../kernel/templateLoader'
import { useTemplateToolStore } from '../../store/templateTool.store'

interface Template {
  id: string
  name: string
  category: string
  data: TemplateData
}

export const useDragAttach = (template: Template) => {
  const { startTemplateDrag } = useTemplateToolStore()

  const handleDragStart = useCallback(
    (ev: React.DragEvent) => {
      ev.dataTransfer.setData('template-id', template.id)
      ev.dataTransfer.effectAllowed = 'copy'
      startTemplateDrag(template.id)
    },
    [template, startTemplateDrag]
  )

  return {
    handleDragStart,
  }
}
