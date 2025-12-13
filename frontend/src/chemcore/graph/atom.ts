export interface Atom {
    id: string;
    // Allowing string to support full periodic table, but validation rules focus on organic subset
    element: string;
    position?: { x: number; y: number; z: number }; // Keeping position for UI compatibility/persistence
}
