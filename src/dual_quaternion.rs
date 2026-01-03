//! Dual quaternions can be used for rigid body movement. They can also be used to represent points and lines.

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
    pub const ZERO:  Self = Self { w: 0.0, i: 0.0, j: 0.0, k: 0.0, ie: 0.0, je: 0.0, ke: 0.0, we: 0.0 };
    pub const ONE:   Self = Self { w:  1.0, ..Self::ZERO };
    pub const DUAL:  Self = Self { we: 1.0, ..Self::ZERO };

    pub const REAL_X: Self = Self { i:  1.0, ..Self::ZERO };
    pub const REAL_Y: Self = Self { j:  1.0, ..Self::ZERO };
    pub const REAL_Z: Self = Self { k:  1.0, ..Self::ZERO };

    pub const DUAL_X: Self = Self { ie:  1.0, ..Self::ZERO };
    pub const DUAL_Y: Self = Self { je:  1.0, ..Self::ZERO };
    pub const DUAL_Z: Self = Self { ke:  1.0, ..Self::ZERO };


    pub fn new(w: Scalar, i: Scalar, j: Scalar, k: Scalar,
               ie: Scalar, je: Scalar, ke: Scalar, we: Scalar) -> Self
    {
        DualQuaternion { w, i, j, k, ie, je, ke, we }
    }

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

    /// Negate the dual part. (= "ideal" or "dual" conjugation)
    pub fn iconj(&self) -> Self
    {
        Self {
            w:   self.w,
            i:   self.i,
            j:   self.j,
            k:   self.k,
            ie: -self.ie,
            je: -self.je,
            ke: -self.ke,
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

    /// The "ideal" or "dual" norm. (= norm of the dual-part-Quaternion)
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
        // So you probably heard of exp(ix) for some real x being some kind of rotational thing.
        // If you multiply some complex number by exp(ix) then it gets rotated.
        //
        // Just imagine we're going 3D now, out of the complex plane.
        // Assume i is just some line through the origin with direction (0,1,0), i.e. directly "into" the complex plane.
        //
        // We can do this rotational thing with all lines, not just through the origin.
        // So if you do exp(ix + Em), for some direction x and moment m, we can create our general rotor and
        // we can even move along the line, so some kind of screwing motion!
        // For screwing we just need to do: exp(ix + Em + Ex).
        //
        // Basically.
        //
        // This algorithm below shows you how to do it exactly.

        let line = line.normalized();
        let angle = 0.5 * angle.rad();
        let distance = 0.5 * distance;

        let dq = DualQuaternion {
            w:  0.0,
            i:  line.i  * angle,
            j:  line.j  * angle,
            k:  line.k  * angle,
            ie: line.ie * angle + distance * line.i,
            je: line.je * angle + distance * line.j,
            ke: line.ke * angle + distance * line.k,
            we: 0.0,
        };

        return dq.exp();
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

    /// Linearily interpolate between `self` and `other`
    pub fn lerp(&self, other: &DualQuaternion, alpha: Scalar) -> DualQuaternion
    {
        (1.0 - alpha) * self + alpha * other
    }

    /// Exponential of a pure dual quaternion.
    /// Will produce wrong results for non-pure dual quaternions.
    /// <div class="warning"> A pure dual quaternion's .we and .w fields are 0.0 <div>
    pub fn exp(&self) -> DualQuaternion
    {
        // You can derive that stuff yourselve if you really want to...
        // It's basically clever grouping and recognizing of terms followed by
        // some power series stuff.
        //
        // You can have a look at this:
        // https://jamessjackson.com/lie_algebra_tutorial/06-closed_form_mat_exp/
        //
        // His $\gamma$ term is my `t` term

        let DualQuaternion { w: _, i, j, k, ie, je, ke, we: _ } = *self;

        let r = (i*i + j*j + k*k).sqrt();

        // exp(0.0) = 1
        // Without this check, it won't work
        if r*r < Scalar::EPSILON {
            return DualQuaternion::ONE
        }

        let t = i*ie + j*je + k*ke;

        let (sin,cos) = r.sin_cos();

        let tr = (cos/(r*r) - sin/(r*r*r))*t;

        DualQuaternion {
            w:   cos,
            i:   (sin/r) * i,
            j:   (sin/r) * j,
            k:   (sin/r) * k,
            ie:  (sin/r) * ie + tr * i,
            je:  (sin/r) * je + tr * j,
            ke:  (sin/r) * ke + tr * k,
            we: -(sin/r) * t
        }
    }

    /// Logarithm of a dual quaternion.
    /// Will produce wrong result when used on unnormalized dual quaternions.
    pub fn log(&self) -> DualQuaternion
    {
        // I took the liberty of taking this algorithm from here, more or less
        // https://jamessjackson.com/lie_algebra_tutorial/06-closed_form_mat_exp/
        //
        // I needed to adjust some minor things though

        let DualQuaternion { w, i, j, k, ie, je, ke, we } = *self;

        let r = (i*i + j*j + k*k).sqrt();
        let t = i*ie + j*je + k*ke;

        let a = (r/w).atan() / r;
        let b = t / (r*r);

        let tr = (w - a)*b - we;

        DualQuaternion {
            w:  0.0,
            i:  a * i,
            j:  a * j,
            k:  a * k,
            ie: a * ie + tr * i,
            je: a * je + tr * j,
            ke: a * ke + tr * k,
            we: 0.0
        }
    }

    /// Raise this dual quaternion to some power.
    /// Will produce incorrect result for unnormalized dual quaternions.
    pub fn powf(&self, f: Scalar) -> DualQuaternion
    {
        ( f * self.log() ).exp()
    }

    /// Screw-lerp this dual quaternion between another dual quaternion.
    /// Only works on normalized dual quaternions.
    pub fn sclerp(&self, other: &DualQuaternion, alpha: Scalar) -> DualQuaternion
    {
        // Took the formula (14) from
        // https://arxiv.org/pdf/2303.13395

        self * (self.conj() * other).powf(alpha)
    }
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
