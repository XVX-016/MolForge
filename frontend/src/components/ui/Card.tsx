import React from 'react';
import clsx from 'clsx';

export type CardProps = {
	header?: React.ReactNode;
	footer?: React.ReactNode;
	children: React.ReactNode;
	className?: string;
};

export default function Card({ header, footer, children, className }: CardProps) {
	return (
		<section
			className={clsx(
				'bg-panel rounded-xl shadow-soft border border-aluminum-DEFAULT',
				className
			)}
		>
			{header && <div className="px-5 py-4 border-b border-aluminum-DEFAULT">{header}</div>}
			<div className="px-5 py-4">{children}</div>
			{footer && <div className="px-5 py-4 border-t border-aluminum-DEFAULT">{footer}</div>}
		</section>
	);
}


