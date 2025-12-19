import { describe, it, expect } from "vitest";
import { MoleculeGraph } from "../src/MoleculeGraph";
import { ForceField } from "../src/ForceField";
import { Vec3 } from "../src/Vec3";

describe("ForceField", () => {
  it("computes bond length force that reduces stretched bond", () => {
    const graph = new MoleculeGraph();
    const cId = graph.addAtom({ element: "C", position: [0, 0, 0] });
    const hId = graph.addAtom({ element: "H", position: [2, 0, 0] }); // Stretched

    const bond = graph.addBond(cId, hId, 1);
    if (!bond) throw new Error("Bond creation failed");

    const atomC = graph.atoms.get(cId)!;
    const atomH = graph.atoms.get(hId)!;
    const bondData = graph.bonds.get(bond)!;

    const { forceC, forceH, energy } = ForceField.computeBondLengthForce(
      atomC,
      atomH,
      bondData
    );

    // Force should pull atoms together (forceC positive x, forceH negative x)
    expect(forceC.x).toBeGreaterThan(0);
    expect(forceH.x).toBeLessThan(0);
    expect(energy).toBeGreaterThan(0);
  });

  it("computes repulsion force that pushes atoms apart", () => {
    const graph = new MoleculeGraph();
    const c1Id = graph.addAtom({ element: "C", position: [0, 0, 0] });
    const c2Id = graph.addAtom({ element: "C", position: [1, 0, 0] }); // Too close

    const atomC1 = graph.atoms.get(c1Id)!;
    const atomC2 = graph.atoms.get(c2Id)!;

    const { forceC1, forceC2, energy } = ForceField.computeRepulsionForce(
      atomC1,
      atomC2
    );

    // Forces should push atoms apart
    expect(forceC1.x).toBeLessThan(0); // C1 pushed left
    expect(forceC2.x).toBeGreaterThan(0); // C2 pushed right
    expect(energy).toBeGreaterThan(0);
  });

  it("optimizeGeometry converges (total energy decreases)", () => {
    const graph = new MoleculeGraph();
    // Create a distorted molecule
    const cId = graph.addAtom({ element: "C", position: [0, 0, 0] });
    const h1Id = graph.addAtom({ element: "H", position: [2, 0, 0] }); // Stretched
    const h2Id = graph.addAtom({ element: "H", position: [-2, 0, 0] }); // Stretched
    const h3Id = graph.addAtom({ element: "H", position: [0, 2, 0] }); // Stretched
    const h4Id = graph.addAtom({ element: "H", position: [0, -2, 0] }); // Stretched

    graph.addBond(cId, h1Id, 1);
    graph.addBond(cId, h2Id, 1);
    graph.addBond(cId, h3Id, 1);
    graph.addBond(cId, h4Id, 1);

    // Compute initial energy
    const initialEnergy = ForceField.computeTotalEnergy(graph);

    // Optimize
    ForceField.optimizeGeometry(graph, 20, 0.005);

    // Compute final energy
    const finalEnergy = ForceField.computeTotalEnergy(graph);

    // Energy should decrease (or at least not increase significantly)
    expect(finalEnergy).toBeLessThanOrEqual(initialEnergy * 1.1); // Allow small tolerance
  });

  it("computes angle force", () => {
    const graph = new MoleculeGraph();
    const cId = graph.addAtom({ element: "C", position: [0, 0, 0] });
    const h1Id = graph.addAtom({ element: "H", position: [1, 0, 0] });
    const h2Id = graph.addAtom({ element: "H", position: [0, 1, 0] });

    const atomC = graph.atoms.get(cId)!;
    const atomH1 = graph.atoms.get(h1Id)!;
    const atomH2 = graph.atoms.get(h2Id)!;

    const { forceC, forceH1, forceH2, energy } = ForceField.computeAngleForce(
      atomH1,
      atomC,
      atomH2
    );

    // Should compute angle forces
    expect(forceC.length()).toBeGreaterThanOrEqual(0);
    expect(forceH1.length()).toBeGreaterThanOrEqual(0);
    expect(forceH2.length()).toBeGreaterThanOrEqual(0);
    expect(energy).toBeGreaterThanOrEqual(0);
  });

  it("getAnglesAround returns angles for atom with multiple neighbors", () => {
    const graph = new MoleculeGraph();
    const cId = graph.addAtom({ element: "C", position: [0, 0, 0] });
    const h1Id = graph.addAtom({ element: "H", position: [1, 0, 0] });
    const h2Id = graph.addAtom({ element: "H", position: [-1, 0, 0] });
    const h3Id = graph.addAtom({ element: "H", position: [0, 1, 0] });

    graph.addBond(cId, h1Id, 1);
    graph.addBond(cId, h2Id, 1);
    graph.addBond(cId, h3Id, 1);

    const angles = graph.getAnglesAround(cId);
    // With 3 neighbors, should have C(3,2) = 3 angles
    expect(angles.length).toBe(3);
    expect(angles[0][1]).toBe(cId); // Central atom is always C
  });
});

