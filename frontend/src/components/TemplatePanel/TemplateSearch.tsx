import React from 'react';

interface TemplateSearchProps {
  search: string;
  setSearch: (value: string) => void;
}

export const TemplateSearch: React.FC<TemplateSearchProps> = ({ search, setSearch }) => {
  return (
    <input
      className="template-search"
      type="text"
      placeholder="Search templates..."
      value={search}
      onChange={e => setSearch(e.target.value)}
    />
  );
};

