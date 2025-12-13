export type BondOrder = 1 | 2 | 3;

export interface Bond {
    id: string;
    from: string;
    to: string;
    order: BondOrder;
}
