import React from 'react';
import clsx from 'clsx';

export type PanelProps = {
	children: React.ReactNode;
	className?: string;
};

export default function Panel({ children, className }: PanelProps) {
	return (
		<div
			className={clsx(
				'bg-panel rounded-xl shadow-soft border border-aluminum-DEFAULT p-4',
				className
			)}
		>
			{children}
		</div>
	);
}


