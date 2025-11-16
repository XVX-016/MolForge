import { useCallback } from 'react';
import { useTemplateToolStore } from '../../store/templateTool.store';

interface Template {
  id: string;
  name: string;
  category: string;
  data: any;
}

export const useDragAttach = (template: Template) => {
  const { startTemplateDrag } = useTemplateToolStore();

  const handleDragStart = useCallback((ev: React.DragEvent) => {
    ev.dataTransfer.setData('template-id', template.id);
    ev.dataTransfer.effectAllowed = 'copy';
    startTemplateDrag(template.id);
  }, [template, startTemplateDrag]);

  return {
    handleDragStart
  };
};

