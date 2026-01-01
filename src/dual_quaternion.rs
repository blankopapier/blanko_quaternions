pub use crate::quaternion::Quaternion;
pub use crate::angle::Angle;

use crate::vector3::Vector3;
use crate::util::Scalar;

#[repr(C)]
#[derive(
    Debug, Clone, Copy, PartialEq, bytemuck::Pod, bytemuck::Zeroable,
    derive_more::Add, derive_more::AddAssign, derive_more::Sub, derive_more::SubAssign,
    derive_more::Neg, derive_more::From
)]
pub struct DualQuaternion
{
    pub w  : Scalar,
    pub i  : Scalar,
    pub j  : Scalar,
    pub k  : Scalar,
    pub ie : Scalar,
    pub je : Scalar,
    pub ke : Scalar,
    pub we : Scalar,
}

impl std::fmt::Display for DualQuaternion
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let components = [
            ("",   self.w),
            ("i",  self.i),
            ("j",  self.j),
            ("k",  self.k),
            ("ie", self.ie),
            ("je", self.je),
            ("ke", self.ke),
            ("we", self.we),
        ];

        for (i,(c,v)) in components.iter().enumerate()
        {
            if v*v <= Scalar::EPSILON {
                continue
            }

            write!(f, "{}{}", c, v);

            if components[i+1..].iter().find(|x| x.1.powi(2) > Scalar::EPSILON).is_some()
            {
                write!(f, " + ");
            }
        }

        write!(f, "")
    }
}

impl DualQuaternion
{
    /// Negate everything except scalar and dual-scalar (clifford conjugation).
    /// Use this when transforming lines in sandwich products
    pub fn conj(&self) -> Self
    {
        Self {
            w:   self.w,
            i:  -self.i,
            j:  -self.j,
            k:  -self.k,
            ie: -self.ie,
            je: -self.je,
            ke: -self.ke,
            we:  self.we,
        }
    }

    /// Same as .conj(), but negates the dual part afterward.
    /// Use this when transforming points in sandwich products
    pub fn nconj(&self) -> Self
    {
        Self {
            w:   self.w,
            i:  -self.i,
            j:  -self.j,
            k:  -self.k,
            ie:  self.ie,
            je:  self.je,
            ke:  self.ke,
            we: -self.we,
        }
    }

    /// Negate the dual part.
    pub fn iconj(&self) -> Self
    {
        Self {
            w:   self.w,
            i:  -self.i,
            j:  -self.j,
            k:  -self.k,
            ie:  self.ie,
            je:  self.je,
            ke:  self.ke,
            we: -self.we,
        }
    }

    /// The norm of the real-part-Quaternion.
    /// For lines, this will the length of the vector3al vector.
    pub fn norm(&self) -> Scalar
    {
        (self.w * self.w +
        self.i * self.i +
        self.j * self.j +
        self.k * self.k).sqrt()
    }

    /// The norm of the dual-part-Quaternion.
    /// For normalized lines, this will the distance from the origin.
    pub fn inorm(&self) -> Scalar
    {
        (self.we * self.we +
        self.ie * self.ie +
        self.je * self.je +
        self.ke * self.ke).sqrt()
    }

    /// Normalize this DualQuaternion by its real-part-Quaternion, i.e. keep rotation normalized.
    /// For lines, this means normalizing the vector3 vector and dividing it's moment by the length
    /// of the original directin vector.
    pub fn normalized(&self) -> Self
    {
        *self * (1.0 / self.norm())
    }

    /// Create a DualQuaternion representing a point in space, i.e. 1+(xi + yj + zk)E.
    pub fn point(pos: &[Scalar]) -> Self
    {
        DualQuaternion { w: 1.0, i: 0.0, j: 0.0, k: 0.0, ie: pos[0], je: pos[1], ke: pos[2], we: 0.0 }
    }

    /// Create a DualQuaternion representing a line in space (= Plücker coordinates)
    ///
    /// Lines can be represented by two points, a point (p) and a vector3 (r) or by a vector3 and a moment (m).
    /// The moment (m) can be calculated with the cross product between a point-vector and a vector3-vector, i.e.
    /// m = p X r
    ///
    /// Let r be the vector (i,j,k) of this DualQuaternion.
    /// Let m be the vector (ie,je,ke) of this DualQuaternion.
    ///
    /// Then this line's nearest point to the origin is P = r X m (cross product of r and m).
    /// Also, r will be normalized.
    pub fn line(pos: &[Scalar], dir: &[Scalar]) -> Self
    {
        let n = 1.0 / (dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]).sqrt();
        let dir = [
            dir[0] * n,
            dir[1] * n,
            dir[2] * n,
        ];

        // normalized cross product between pos and dir
        let moment = [
            n * (pos[1] * dir[2] - dir[1] * pos[2]),
            n * (pos[2] * dir[0] - dir[2] * pos[0]),
            n * (pos[0] * dir[1] - dir[0] * pos[1]),
        ];

        DualQuaternion { w: 0.0, i: dir[0], j: dir[1], k: dir[2], ie: moment[0], je: moment[1], ke: moment[2], we: 0.0 }
    }

    /// Create a screw around a line.
    /// The screw will rotate `angle` and travel `distance` units along the line.
    pub fn screw(line: &DualQuaternion, angle: Angle, distance: Scalar) -> Self
    {
        // So basically, we use Geometric Algebra for this and just call it DualQuaternion afterward.
        // The key idea is to convert position and vector3 to moment and vector3, aka Plücker Coordinates for a line.
        // We then 'normalize' this line and create a rotor and translator from it.
        // Then compose the two DualQuaternions and return it.

        let (sin,cos) = (angle*0.5).sin_cos();
        let t = distance * 0.5;

        let rotor = DualQuaternion {
            w:  cos,
            i:  sin * line.i,
            j:  sin * line.j,
            k:  sin * line.k,
            ie: sin * line.ie,
            je: sin * line.je,
            ke: sin * line.ke,
            we: 0.0
        };

        let translator = DualQuaternion {
            w:  1.0,
            i:  0.0,
            j:  0.0,
            k:  0.0,
            ie: t * line.i,
            je: t * line.j,
            ke: t * line.k,
            we: 0.0
        };

        translator * rotor
    }

    /// Basically a Quaternion
    pub fn rotor(angle: Angle, axis: &[Scalar]) -> Self
    {
        let axis = Vector3 {
            x: axis[0],
            y: axis[1],
            z: axis[2]
        };

        let (sin,cos) = (angle*0.5).sin_cos();
        let dir = axis.normalize();

        DualQuaternion {
            w:  cos,
            i:  sin * dir.x,
            j:  sin * dir.y,
            k:  sin * dir.z,
            ie: 0.0,
            je: 0.0,
            ke: 0.0,
            we: 0.0
        }
    }

    /// Translational DualQuaternion.
    pub fn translator(translation: &[Scalar]) -> Self
    {
        DualQuaternion {
            w:  1.0,
            i:  0.0,
            j:  0.0,
            k:  0.0,
            ie: 0.5 * translation[0],
            je: 0.5 * translation[1],
            ke: 0.5 * translation[2],
            we: 0.0
        }
    }

    /// Transform a 3D-vector as point.
    /// This means that the vector will be screwed around a line.
    pub fn transform_point(&self, point: &[Scalar]) -> [Scalar; 3]
    {
        let point = Vector3 {
            x: point[0],
            y: point[1],
            z: point[2]
        };

        // Taken from
        // https://rigidgeometricalgebra.org/wiki/index.php?title=Motor

        let v = Vector3 { x: self.i,  y: self.j,  z: self.k  };
        let m = Vector3 { x: self.ie, y: self.je, z: self.ke };
        let vw = self.w;
        let mw = self.we;

        let p: Vector3 = point.into();
        let a = v.cross(&p) + m;

        (p + 2.0 * (vw*a + v.cross(&a) - mw*v)).into()
    }

    /// Transform a 3D-vector as vector3.
    /// This means that the vector will be rotated around the origin, not
    /// around a line. Neither will it be translated along a line.
    pub fn transform_vector3(&self, vector3: &[Scalar]) -> [Scalar; 3]
    {
        let vector3 = Vector3 {
            x: vector3[0],
            y: vector3[1],
            z: vector3[2]
        };

        // Taken from
        // https://rigidgeometricalgebra.org/wiki/index.php?title=Motor

        let v = Vector3 { x: self.i,  y: self.j,  z: self.k  };
        let vw = self.w;

        let a = v.cross(&vector3);

        (vector3 + 2.0 * (vw*a + v.cross(&a))).into()
    }

    /// Transform a line by this DualQuaternion, i.e. screw some line around another line
    pub fn transform_line(&self, line: &DualQuaternion) -> DualQuaternion
    {
        // Taken from
        // https://rigidgeometricalgebra.org/wiki/index.php?title=Motor
        //
        // This is the same as transforming vector3 and position of the line
        // individually and converting it to a DualQuaternion again.

        let lv = Vector3 { x: line.i, y: line.j, z: line.k };
        let lm = Vector3 { x: line.ie, y: line.je, z: line.ke };

        let v = Vector3 { x: self.i,  y: self.j,  z: self.k  };
        let m = Vector3 { x: self.ie, y: self.je, z: self.ke };
        let vw = self.w;
        let mw = self.we;

        let a = v.cross(&lv);
        let b = v.cross(&lm);
        let c = m.cross(&lv);
        let d = b + c; // Not

        let lv = lv + 2.0 * (vw*a + v.cross(&a));
        let lm = lm + 2.0 * (mw*a + vw*d + v.cross(&d) + m.cross(&a));

        DualQuaternion {
            w: 0.0,
            i: lv.x, j: lv.y, k: lv.z,
            ie: lm.x, je: lm.y, ke: lm.z,
            we: 0.0
        }
    }

    // TODO: Pow, Log, Exp
}

auto_ops::impl_op_ex!(* |lhs: &DualQuaternion, rhs: &DualQuaternion| -> DualQuaternion {
    let (lhs_real, lhs_dual) = (
        Quaternion { w: lhs.w, i: lhs.i, j: lhs.j, k: lhs.k },
        Quaternion { w: lhs.we, i: lhs.ie, j: lhs.je, k: lhs.ke }
    );

    let (rhs_real, rhs_dual) = (
        Quaternion { w: rhs.w, i: rhs.i, j: rhs.j, k: rhs.k },
        Quaternion { w: rhs.we, i: rhs.ie, j: rhs.je, k: rhs.ke }
    );

    let (Quaternion { w, i, j, k }, Quaternion { w: we, i: ie, j: je, k: ke }) = (
        lhs_real * rhs_real,
        lhs_real * rhs_dual + lhs_dual * rhs_real
    );

    DualQuaternion { w, i, j, k, ie, je, ke, we }
});

auto_ops::impl_op_ex_commutative!(* |lhs: &DualQuaternion, rhs: &Scalar| -> DualQuaternion {
    DualQuaternion {
        w:  rhs * lhs.w,
        i:  rhs * lhs.i,
        j:  rhs * lhs.j,
        k:  rhs * lhs.k,
        ie: rhs * lhs.ie,
        je: rhs * lhs.je,
        ke: rhs * lhs.ke,
        we: rhs * lhs.we
    }
});

auto_ops::impl_op_ex!(*= |lhs: &mut DualQuaternion, rhs: &Scalar| {
    lhs.w  = lhs.w  * rhs;
    lhs.i  = lhs.i  * rhs;
    lhs.j  = lhs.j  * rhs;
    lhs.k  = lhs.k  * rhs;
    lhs.ie = lhs.ie * rhs;
    lhs.je = lhs.je * rhs;
    lhs.ke = lhs.ke * rhs;
    lhs.we = lhs.we * rhs;
});
