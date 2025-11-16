import React, { useState, useMemo } from 'react';
import { getTemplates } from '../../templates/templates.index';
import { TemplateCategory } from './TemplateCategory';
import { TemplateSearch } from './TemplateSearch';
import './TemplatePanel.css';

export const TemplatePanel: React.FC = () => {
  const templates = getTemplates();
  const [search, setSearch] = useState('');

  const filtered = useMemo(() => {
    if (!search.trim()) return templates;
    return templates.filter(t =>
      t.name.toLowerCase().includes(search.toLowerCase()) ||
      t.category.toLowerCase().includes(search.toLowerCase())
    );
  }, [search, templates]);

  const categories = useMemo(() => {
    const map = new Map<string, typeof templates>();
    for (const t of filtered) {
      if (!map.has(t.category)) map.set(t.category, []);
      map.get(t.category)!.push(t);
    }
    return Array.from(map.entries());
  }, [filtered]);

  return (
    <div className="template-panel">
      <TemplateSearch search={search} setSearch={setSearch} />
      {categories.map(([cat, items]) => (
        <TemplateCategory key={cat} category={cat} items={items} />
      ))}
    </div>
  );
};
