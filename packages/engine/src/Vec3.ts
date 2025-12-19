/**
 * 3D Vector utility class for molecular mechanics calculations
 */
export class Vec3 {
  constructor(
    public x: number = 0,
    public y: number = 0,
    public z: number = 0
  ) {}

  /**
   * Create from array [x, y, z]
   */
  static fromArray(arr: [number, number, number]): Vec3 {
    return new Vec3(arr[0], arr[1], arr[2])
  }

  /**
   * Create from position tuple
   */
  static fromPosition(pos: [number, number, number]): Vec3 {
    return new Vec3(pos[0], pos[1], pos[2])
  }

  /**
   * Convert to array [x, y, z]
   */
  toArray(): [number, number, number] {
    return [this.x, this.y, this.z]
  }

  /**
   * Add another vector
   */
  add(other: Vec3): Vec3 {
    return new Vec3(this.x + other.x, this.y + other.y, this.z + other.z)
  }

  /**
   * Subtract another vector
   */
  sub(other: Vec3): Vec3 {
    return new Vec3(this.x - other.x, this.y - other.y, this.z - other.z)
  }

  /**
   * Multiply by scalar
   */
  mul(scalar: number): Vec3 {
    return new Vec3(this.x * scalar, this.y * scalar, this.z * scalar)
  }

  /**
   * Divide by scalar
   */
  div(scalar: number): Vec3 {
    return new Vec3(this.x / scalar, this.y / scalar, this.z / scalar)
  }

  /**
   * Dot product
   */
  dot(other: Vec3): number {
    return this.x * other.x + this.y * other.y + this.z * other.z
  }

  /**
   * Cross product
   */
  cross(other: Vec3): Vec3 {
    return new Vec3(
      this.y * other.z - this.z * other.y,
      this.z * other.x - this.x * other.z,
      this.x * other.y - this.y * other.x
    )
  }

  /**
   * Vector length (magnitude)
   */
  length(): number {
    return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z)
  }

  /**
   * Squared length (for performance, avoids sqrt)
   */
  lengthSq(): number {
    return this.x * this.x + this.y * this.y + this.z * this.z
  }

  /**
   * Normalize to unit vector
   */
  normalize(): Vec3 {
    const len = this.length()
    if (len === 0) return new Vec3(0, 0, 0)
    return this.div(len)
  }

  /**
   * Angle between two vectors (in radians)
   */
  static angleBetween(a: Vec3, b: Vec3): number {
    const dot = a.dot(b)
    const lenA = a.length()
    const lenB = b.length()
    if (lenA === 0 || lenB === 0) return 0
    return Math.acos(Math.max(-1, Math.min(1, dot / (lenA * lenB))))
  }

  /**
   * Distance between two vectors
   */
  static distance(a: Vec3, b: Vec3): number {
    return a.sub(b).length()
  }

  /**
   * Copy vector
   */
  clone(): Vec3 {
    return new Vec3(this.x, this.y, this.z)
  }
}

