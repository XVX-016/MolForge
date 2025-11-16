import React, { useState } from 'react';
import { TemplateItem } from './TemplateItem';

interface Template {
  id: string;
  name: string;
  category: string;
  data: any;
}

interface TemplateCategoryProps {
  category: string;
  items: Template[];
}

export const TemplateCategory: React.FC<TemplateCategoryProps> = ({ category, items }) => {
  const [open, setOpen] = useState(true);

  return (
    <div className="template-category">
      <div className="template-category-header" onClick={() => setOpen(o => !o)}>
        <span>{category}</span>
        <span className="count">{items.length}</span>
      </div>

      {open && (
        <div className="template-category-grid">
          {items.map(t => (
            <TemplateItem key={t.id} template={t} />
          ))}
        </div>
      )}
    </div>
  );
};
