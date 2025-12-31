use crate::quaternion::Quaternion;
pub use crate::angle::Angle;
pub use crate::point::Point;
pub use crate::direction::Direction;

#[repr(C)]
#[derive(
    Debug, Clone, Copy, PartialEq, bytemuck::Pod, bytemuck::Zeroable,
    derive_more::Add, derive_more::AddAssign, derive_more::Sub, derive_more::SubAssign,
    derive_more::Neg, derive_more::From
)]
pub struct DualQuaternion
{
    pub w  : f32,
    pub i  : f32,
    pub j  : f32,
    pub k  : f32,
    pub ie : f32,
    pub je : f32,
    pub ke : f32,
    pub we : f32,
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
            if v*v <= std::f32::EPSILON {
                continue
            }

            write!(f, "{}{}", c, v);

            if components[i+1..].iter().find(|x| x.1.powi(2) > std::f32::EPSILON).is_some()
            {
                write!(f, " + ");
            }
        }

        write!(f, "")
    }
}

impl DualQuaternion
{
    /// Negate everything except scalar and dual-scalar (clifford conjugation)
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

    /// "Sandwich"-Conjugate, use this with sandwich-product.
    /// Just like DualQuaternion::conj(), but negate the dual part afterward
    pub fn sconj(&self) -> Self
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

    /// The norm of the real-part-Quaternion
    pub fn norm(&self) -> f32
    {
        (self.w * self.w +
        self.i * self.i +
        self.j * self.j +
        self.k * self.k).sqrt()
    }

    /// The norm of the dual-part-Quaternion
    pub fn inorm(&self) -> f32
    {
        (self.we * self.we +
        self.ie * self.ie +
        self.je * self.je +
        self.ke * self.ke).sqrt()
    }

    /// Normalize this DualQuaternion by its real-part-Quaternion, i.e. keep rotation normalized
    pub fn normalized(&self) -> Self
    {
        *self * (1.0 / self.norm())
    }

    /// Screw around a line through `pos` with direction `dir`.
    /// The screw will rotate `angle` and travel `distance` units along the line.
    pub fn from_line(pos: &Point, dir: &Direction, angle: Angle, distance: f32) -> Self
    {
        // So basically, we use Geometric Algebra for this and just call it DualQuaternion afterward.
        // The key idea is to convert position and direction to moment and direction, aka Plücker Coordinates for a line.
        // We then 'normalize' this line and create a rotor and translator from it.
        // Then compose the two DualQuaternions and return it.

        let normalizer = 1.0 / dir.norm();
        let (sin,cos) = (angle*0.5).sin_cos();
        let t = distance * 0.5;

        let pos: Direction = pos.into();
        let moment = pos.cross(dir);

        let dir = dir * normalizer;
        let moment = moment * normalizer;

        let rotor = DualQuaternion {
            w:  cos,
            i:  sin * dir.x,
            j:  sin * dir.y,
            k:  sin * dir.z,
            ie: sin * moment.x,
            je: sin * moment.y,
            ke: sin * moment.z,
            we: 0.0
        };

        let translator = DualQuaternion {
            w:  1.0,
            i:  0.0,
            j:  0.0,
            k:  0.0,
            ie: t * dir.x,
            je: t * dir.y,
            ke: t * dir.z,
            we: 0.0
        };

        translator * rotor
    }

    /// Basically a Quaternion
    pub fn from_angle_axis(angle: Angle, axis: &Direction) -> Self
    {
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
    pub fn from_translation(translation: &Direction) -> Self
    {
        DualQuaternion {
            w:  1.0,
            i:  0.0,
            j:  0.0,
            k:  0.0,
            ie: 0.5 * translation.x,
            je: 0.5 * translation.y,
            ke: 0.5 * translation.z,
            we: 0.0
        }
    }

    pub fn transform_point(&self, point: &Point) -> Point
    {
        // Taken from
        // https://rigidgeometricalgebra.org/wiki/index.php?title=Motor

        let v = Direction { x: self.i,  y: self.j,  z: self.k  };
        let m = Direction { x: self.ie, y: self.je, z: self.ke };
        let vw = self.w;
        let mw = self.we;

        let p: Direction = point.into();
        let a = v.cross(&p) + m;

        (p + 2.0 * (vw*a + v.cross(&a) - mw*v)).into()
    }

    pub fn transform_direction(&self, direction: &Direction) -> Direction
    {
        // Taken from
        // https://rigidgeometricalgebra.org/wiki/index.php?title=Motor

        let v = Direction { x: self.i,  y: self.j,  z: self.k  };
        let vw = self.w;

        let a = v.cross(&direction);

        *direction + 2.0 * (vw*a + v.cross(&a))
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
auto_ops::impl_op_ex!(* |lhs: &DualQuaternion, rhs: &Quaternion| -> DualQuaternion {
    let (lhs_real, lhs_dual) = (
        Quaternion { w: lhs.w, i: lhs.i, j: lhs.j, k: lhs.k },
        Quaternion { w: lhs.we, i: lhs.ie, j: lhs.je, k: lhs.ke }
    );

    let (Quaternion { w, i, j, k }, Quaternion { w: we, i: ie, j: je, k: ke }) = (
        lhs_real * rhs,
        lhs_dual * rhs
    );

    DualQuaternion { w, i, j, k, ie, je, ke, we }
});
auto_ops::impl_op_ex!(* |lhs: &Quaternion, rhs: &DualQuaternion| -> DualQuaternion {
    let (rhs_real, rhs_dual) = (
        Quaternion { w: rhs.w, i: rhs.i, j: rhs.j, k: rhs.k },
        Quaternion { w: rhs.we, i: rhs.ie, j: rhs.je, k: rhs.ke }
    );

    let (Quaternion { w, i, j, k }, Quaternion { w: we, i: ie, j: je, k: ke }) = (
        lhs * rhs_real,
        lhs * rhs_dual
    );

    DualQuaternion { w, i, j, k, ie, je, ke, we }
});
auto_ops::impl_op_ex_commutative!(* |lhs: &DualQuaternion, rhs: &f32| -> DualQuaternion {
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
auto_ops::impl_op_ex!(*= |lhs: &mut DualQuaternion, rhs: &f32| {
    lhs.w  = lhs.w  * rhs;
    lhs.i  = lhs.i  * rhs;
    lhs.j  = lhs.j  * rhs;
    lhs.k  = lhs.k  * rhs;
    lhs.ie = lhs.ie * rhs;
    lhs.je = lhs.je * rhs;
    lhs.ke = lhs.ke * rhs;
    lhs.we = lhs.we * rhs;
});
