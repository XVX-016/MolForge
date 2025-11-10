import { describe, it, expect } from 'vitest';
import React from 'react';
import { renderToString } from 'react-dom/server';
import Button from '../../ui/Button';

describe('Button', () => {
	it('renders primary by default', () => {
		const html = renderToString(<Button>Click</Button>);
		expect(html).toContain('Click');
	});

	it('renders secondary variant', () => {
		const html = renderToString(<Button variant="secondary">Go</Button>);
		expect(html).toContain('Go');
	});
});


