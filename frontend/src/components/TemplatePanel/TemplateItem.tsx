import React from 'react';
import { placeTemplate } from '../../kernel/templateLoader';
import { useDragAttach } from './useDragAttach';

interface Template {
  id: string;
  name: string;
  category: string;
  data: any;
}

interface TemplateItemProps {
  template: Template;
}

export const TemplateItem: React.FC<TemplateItemProps> = ({ template }) => {
  const { handleDragStart } = useDragAttach(template);

  const handleClick = () => {
    placeTemplate(template.data, { x: 0, y: 0, z: 0 });
  };

  return (
    <div
      className="template-item"
      draggable
      onDragStart={handleDragStart}
      onClick={handleClick}
    >
      <div className="template-icon" />
      <div className="template-name">{template.name}</div>
    </div>
  );
};
