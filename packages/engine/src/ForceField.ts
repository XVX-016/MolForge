import { MoleculeGraph } from "./MoleculeGraph";
import { Atom, Bond, getBondLength } from "./types";
import { Vec3 } from "./Vec3";

/**
 * Force field parameters
 */
const FORCE_CONSTANTS = {
  bond: {
    "C-C": 350.0, // kcal/mol/Å²
    "C-H": 340.0,
    "C-O": 360.0,
    "C-N": 320.0,
    "O-H": 450.0,
    "N-H": 340.0,
    default: 300.0,
  },
  angle: {
    default: 40.0, // kcal/mol/rad²
  },
};

const LENNARD_JONES = {
  epsilon: 0.1, // kcal/mol
  sigma: 3.0, // Å
};

/**
 * Simplified molecular mechanics force field
 */
export class ForceField {
  /**
   * Compute bond stretching force
   * E_stretch = k_bond * (d - d0)^2
   */
  static computeBondLengthForce(
    atomA: Atom,
    atomB: Atom,
    bond: Bond
  ): { forceA: Vec3; forceB: Vec3; energy: number } {
    const posA = Vec3.fromPosition(atomA.position);
    const posB = Vec3.fromPosition(atomB.position);
    const vec = posB.sub(posA);
    const distance = vec.length();
    const desiredLength = getBondLength(atomA.element, atomB.element);

    // Force constant based on bond type
    const bondKey = `${atomA.element}-${atomB.element}`;
    const k_bond =
      FORCE_CONSTANTS.bond[bondKey as keyof typeof FORCE_CONSTANTS.bond] ||
      FORCE_CONSTANTS.bond.default;

    // Energy: E = k * (d - d0)^2
    const energy = k_bond * Math.pow(distance - desiredLength, 2);

    // Force magnitude: F = -dE/dr = -2k(d - d0)
    const forceMagnitude = -2 * k_bond * (distance - desiredLength);

    if (distance > 0.001) {
      const direction = vec.normalize();
      const force = direction.mul(forceMagnitude);

      // Forces are equal and opposite
      return {
        forceA: force,
        forceB: force.mul(-1),
        energy,
      };
    }

    return {
      forceA: new Vec3(0, 0, 0),
      forceB: new Vec3(0, 0, 0),
      energy: 0,
    };
  }

  /**
   * Compute angle bending force
   * E_angle = k_angle * (theta - theta0)^2
   */
  static computeAngleForce(
    atomA: Atom,
    atomB: Atom,
    atomC: Atom
  ): { forceA: Vec3; forceB: Vec3; forceC: Vec3; energy: number } {
    const posA = Vec3.fromPosition(atomA.position);
    const posB = Vec3.fromPosition(atomB.position);
    const posC = Vec3.fromPosition(atomC.position);

    const vecBA = posA.sub(posB);
    const vecBC = posC.sub(posB);

    const lenBA = vecBA.length();
    const lenBC = vecBC.length();

    if (lenBA < 0.001 || lenBC < 0.001) {
      return {
        forceA: new Vec3(0, 0, 0),
        forceB: new Vec3(0, 0, 0),
        forceC: new Vec3(0, 0, 0),
        energy: 0,
      };
    }

    // Compute angle
    const cosTheta = vecBA.dot(vecBC) / (lenBA * lenBC);
    const theta = Math.acos(Math.max(-1, Math.min(1, cosTheta)));

    // Desired angle (tetrahedral for most organic molecules)
    const theta0 = (109.5 * Math.PI) / 180; // ~109.5 degrees in radians

    const k_angle = FORCE_CONSTANTS.angle.default;

    // Energy: E = k * (theta - theta0)^2
    const energy = k_angle * Math.pow(theta - theta0, 2);

    // Force magnitude: F = -dE/dtheta = -2k(theta - theta0)
    const forceMagnitude = -2 * k_angle * (theta - theta0);

    // Simplified force distribution (perpendicular to bond vectors)
    const normBA = vecBA.normalize();
    const normBC = vecBC.normalize();
    const perpA = normBA.sub(normBC.mul(normBA.dot(normBC))).normalize();
    const perpC = normBC.sub(normBA.mul(normBC.dot(normBA))).normalize();

    return {
      forceA: perpA.mul(forceMagnitude / lenBA),
      forceB: perpA.mul(-forceMagnitude / lenBA).add(perpC.mul(-forceMagnitude / lenBC)),
      forceC: perpC.mul(forceMagnitude / lenBC),
      energy,
    };
  }

  /**
   * Compute Lennard-Jones repulsion force
   * E_lj = 4 * epsilon * ( (sigma/r)^12 - (sigma/r)^6 )
   */
  static computeRepulsionForce(
    atomA: Atom,
    atomB: Atom
  ): { forceA: Vec3; forceB: Vec3; energy: number } {
    const posA = Vec3.fromPosition(atomA.position);
    const posB = Vec3.fromPosition(atomB.position);
    const vec = posB.sub(posA);
    const distance = vec.length();

    if (distance < 0.001) {
      return {
        forceA: new Vec3(0, 0, 0),
        forceB: new Vec3(0, 0, 0),
        energy: 0,
      };
    }

    const epsilon = LENNARD_JONES.epsilon;
    const sigma = LENNARD_JONES.sigma;

    const r6 = Math.pow(sigma / distance, 6);
    const r12 = r6 * r6;

    // Energy: E = 4 * epsilon * (r12 - r6)
    const energy = 4 * epsilon * (r12 - r6);

    // Force magnitude: F = -dE/dr = (24 * epsilon / r) * (2*r12 - r6)
    const forceMagnitude = (24 * epsilon / distance) * (2 * r12 - r6);

    const direction = vec.normalize();
    const force = direction.mul(forceMagnitude);

    return {
      forceA: force.mul(-1),
      forceB: force,
      energy,
    };
  }

  /**
   * Compute total energy of molecule
   */
  static computeTotalEnergy(molecule: MoleculeGraph): number {
    let totalEnergy = 0;

    // Bond stretching energy
    molecule.bonds.forEach((bond) => {
      const atomA = molecule.atoms.get(bond.a1);
      const atomB = molecule.atoms.get(bond.a2);
      if (atomA && atomB) {
        const { energy } = this.computeBondLengthForce(atomA, atomB, bond);
        totalEnergy += energy;
      }
    });

    // Angle bending energy
    const angles = molecule.getAllAngles();
    angles.forEach(([idA, idB, idC]) => {
      const atomA = molecule.atoms.get(idA);
      const atomB = molecule.atoms.get(idB);
      const atomC = molecule.atoms.get(idC);
      if (atomA && atomB && atomC) {
        const { energy } = this.computeAngleForce(atomA, atomB, atomC);
        totalEnergy += energy;
      }
    });

    // Non-bonded repulsion (only for non-bonded pairs)
    const atoms = Array.from(molecule.atoms.values());
    for (let i = 0; i < atoms.length; i++) {
      for (let j = i + 1; j < atoms.length; j++) {
        const atomA = atoms[i];
        const atomB = atoms[j];

        // Check if atoms are bonded
        const isBonded = Array.from(molecule.bonds.values()).some(
          (b) =>
            (b.a1 === atomA.id && b.a2 === atomB.id) ||
            (b.a1 === atomB.id && b.a2 === atomA.id)
        );

        if (!isBonded) {
          const { energy } = this.computeRepulsionForce(atomA, atomB);
          totalEnergy += energy;
        }
      }
    }

    return totalEnergy;
  }

  /**
   * Optimize geometry using gradient descent
   */
  static optimizeGeometry(
    molecule: MoleculeGraph,
    iterations: number = 20,
    step: number = 0.005
  ): void {
    const atoms = Array.from(molecule.atoms.values());

    for (let iter = 0; iter < iterations; iter++) {
      // Initialize forces for all atoms
      const forces = new Map<string, Vec3>();
      atoms.forEach((atom) => {
        forces.set(atom.id, new Vec3(0, 0, 0));
      });

      // Compute bond stretching forces
      molecule.bonds.forEach((bond) => {
        const atomA = molecule.atoms.get(bond.a1);
        const atomB = molecule.atoms.get(bond.a2);
        if (atomA && atomB) {
          const { forceA, forceB } = this.computeBondLengthForce(
            atomA,
            atomB,
            bond
          );
          const fA = forces.get(atomA.id)!;
          const fB = forces.get(atomB.id)!;
          forces.set(atomA.id, fA.add(forceA));
          forces.set(atomB.id, fB.add(forceB));
        }
      });

      // Compute angle bending forces
      const angles = molecule.getAllAngles();
      angles.forEach(([idA, idB, idC]) => {
        const atomA = molecule.atoms.get(idA);
        const atomB = molecule.atoms.get(idB);
        const atomC = molecule.atoms.get(idC);
        if (atomA && atomB && atomC) {
          const { forceA, forceB, forceC } = this.computeAngleForce(
            atomA,
            atomB,
            atomC
          );
          const fA = forces.get(atomA.id)!;
          const fB = forces.get(atomB.id)!;
          const fC = forces.get(atomC.id)!;
          forces.set(atomA.id, fA.add(forceA));
          forces.set(atomB.id, fB.add(forceB));
          forces.set(atomC.id, fC.add(forceC));
        }
      });

      // Compute non-bonded repulsion forces
      for (let i = 0; i < atoms.length; i++) {
        for (let j = i + 1; j < atoms.length; j++) {
          const atomA = atoms[i];
          const atomB = atoms[j];

          // Check if atoms are bonded
          const isBonded = Array.from(molecule.bonds.values()).some(
            (b) =>
              (b.a1 === atomA.id && b.a2 === atomB.id) ||
              (b.a1 === atomB.id && b.a2 === atomA.id)
          );

          if (!isBonded) {
            const { forceA, forceB } = this.computeRepulsionForce(atomA, atomB);
            const fA = forces.get(atomA.id)!;
            const fB = forces.get(atomB.id)!;
            forces.set(atomA.id, fA.add(forceA));
            forces.set(atomB.id, fB.add(forceB));
          }
        }
      }

      // Update positions using gradient descent
      atoms.forEach((atom) => {
        const force = forces.get(atom.id)!;
        const displacement = force.mul(step);
        const newPos = Vec3.fromPosition(atom.position).add(displacement);
        molecule.atoms.set(atom.id, {
          ...atom,
          position: newPos.toArray(),
        });
      });
    }
  }
}

