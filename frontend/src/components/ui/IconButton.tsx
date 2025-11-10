import React from 'react';
import clsx from 'clsx';

export type IconButtonProps = React.ButtonHTMLAttributes<HTMLButtonElement> & {
	'aria-label': string;
};

export default function IconButton({ className, children, ...rest }: IconButtonProps) {
	return (
		<button
			type="button"
			{...rest}
			className={clsx(
				'inline-flex h-9 w-9 items-center justify-center rounded-lg border border-aluminum-DEFAULT bg-panel text-text-secondary hover:text-text-primary hover:bg-aluminum-light focus:outline-none focus:ring-2 focus:ring-accent-blue',
				className
			)}
		>
			{children}
		</button>
	);
}


