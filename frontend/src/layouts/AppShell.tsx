import React from 'react';
import Navbar from '../components/Navbar';

type AppShellProps = {
	children: React.ReactNode;
	noPadding?: boolean;
};

export default function AppShell({ children, noPadding = false }: AppShellProps) {
	return (
		<div className="h-screen bg-white text-black flex flex-col overflow-hidden">
			<Navbar />
			<main className={`flex-1 flex flex-col overflow-y-auto overflow-x-hidden ${noPadding ? '' : 'p-4 sm:p-6 lg:p-8'}`}>
				{children}
			</main>
		</div>
	);
}
