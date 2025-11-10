import React from 'react';
import clsx from 'clsx';

export type ButtonProps = React.ButtonHTMLAttributes<HTMLButtonElement> & {
	variant?: 'primary' | 'secondary';
};

export default function Button({
	variant = 'primary',
	className,
	disabled,
	children,
	...rest
}: ButtonProps) {
	const classes =
		variant === 'primary'
			? 'bg-accent-blue text-white hover:brightness-95 disabled:opacity-50'
			: 'bg-aluminum-light text-text-primary hover:bg-aluminum-dark/30 disabled:opacity-50';

	return (
		<button
			{...rest}
			disabled={disabled}
			className={clsx(
				'inline-flex items-center justify-center rounded-lg px-4 py-2 font-medium transition-colors focus:outline-none focus:ring-2 focus:ring-accent-blue',
				classes,
				className
			)}
		>
			{children}
		</button>
	);
}


