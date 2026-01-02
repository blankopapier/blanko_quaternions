// Quaternions


pub use crate::angle::Angle;
use crate::vector3::Vector3;
use crate::util::Scalar;

#[repr(C)]
#[derive(
    Debug, Clone, Copy, PartialEq, bytemuck::Pod, bytemuck::Zeroable,
    derive_more::Add, derive_more::AddAssign, derive_more::Sub, derive_more::SubAssign,
    derive_more::Neg, derive_more::From
)]
pub struct Quaternion
{
    pub w: Scalar,
    pub i: Scalar,
    pub j: Scalar,
    pub k: Scalar,
}

impl std::fmt::Display for Quaternion
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.w.powi(2) > Scalar::EPSILON {
            write!(f, "{}",  self.w);

            if self.i.powi(2) > Scalar::EPSILON ||
                self.j.powi(2) > Scalar::EPSILON ||
                self.k.powi(2) > Scalar::EPSILON
            {
                write!(f, " + ");
            }
        }

        if self.i.powi(2) > Scalar::EPSILON {
            write!(f, "{}i", self.i);

            if self.j.powi(2) > Scalar::EPSILON ||
                self.k.powi(2) > Scalar::EPSILON
            {
                write!(f, " + ");
            }
        }

        if self.j.powi(2) > Scalar::EPSILON {
            write!(f, "{}j", self.j);

            if self.k.powi(2) > Scalar::EPSILON
            {
                write!(f, " + ");
            }
        }

        if self.k.powi(2) > Scalar::EPSILON {
            write!(f, "{}k", self.k);
        }

        write!(f, "")
    }
}

impl Quaternion
{
    pub fn conj(&self) -> Self { Self { w: self.w, i: -self.i, j: -self.j, k: -self.k } }
    pub fn norm(&self) -> Scalar { (self.w*self.w + self.i*self.i + self.j*self.j + self.k*self.k).sqrt() }
    pub fn normalized(&self) -> Self { *self * (1.0 / self.norm()) }

    /// Get this Quaternion's angle, i.e. if q=r(cos(x)+sin(x)*v), where r is the norm and v is a normalized axis, get x
    pub fn angle(&self) -> Angle
    {
        let n = self.normalized();
        Angle::radians(n.w.acos())
    }

    /// Create a Quaternion representing a point in space, i.e. xi + yj + zk.
    pub fn point(pos: &[Scalar]) -> Self
    {
        Quaternion { w: 0.0, i: pos[0], j: pos[1], k: pos[2] }
    }

    pub fn rotor(angle: Angle, x: Scalar, y: Scalar, z: Scalar) -> Self
    {
        let mut q = Self { w: 0.0, i: x, j: y, k: z }.normalized();
        let (sin,cos) = (angle*0.5).sin_cos();

        q *= sin;
        q.w = cos;

        q
    }

    pub fn transform_point(&self, point: &[Scalar]) -> [Scalar;3]
    {
        self.transform_direction(point).into()
    }

    pub fn transform_direction(&self, direction: &[Scalar]) -> [Scalar;3]
    {
        // Taken from
        // https://rigidgeometricalgebra.org/wiki/index.php?title=Motor

        let direction = Vector3 {
            x: direction[0],
            y: direction[1],
            z: direction[2]
        };

        let v = Vector3 { x: self.i,  y: self.j,  z: self.k  };
        let vw = self.w;

        let a = v.cross(&direction);

        (direction + 2.0 * (vw*a + v.cross(&a))).into()
    }

    /// Linearily interpolate between this and `other`
    pub fn lerp(&self, other: Quaternion, alpha: Scalar) -> Quaternion
    {
        (1.0 - alpha) * self + alpha * other
    }

    /// Spherically interpolate between `self` and `other`
    pub fn slerp(&self, other: Quaternion, alpha: Scalar) -> Quaternion
    {
        // https://en.wikipedia.org/wiki/Slerp
        // Extended to interpolate length as well

        let (r1,r2) = (self.norm(), other.norm());
        let (q1,q2) = (self.normalized(), other.normalized());

        let q = ( q2*q1.conj() ).powf(alpha) * q1;

        ( (1.0 - alpha) * r1 + alpha * r2 ) * q
    }

    /// Raise this Quaternion to a (real) power
    pub fn powf(&self, f: Scalar) -> Quaternion
    {
        let r = self.norm();

        // Check if this Quaternion is just a scalar
        if (self.w - r).powi(2) <= Scalar::EPSILON
        {
            Quaternion { w: r.powf(f), i: 0.0, j: 0.0, k: 0.0 }
        }
        else
        {
            let n    = self * (1.0/r);
            let axis = Vector3 { x: n.i, y: n.j, z: n.k }.normalize();

            let angle     = n.w.acos() * f;
            let (sin,cos) = angle.sin_cos();

            let r = r.powf(f);

            Quaternion {
                w: r * cos,
                i: r * sin * axis.x,
                j: r * sin * axis.y,
                k: r * sin * axis.z,
            }

        }

    }

    /// Exponential of a quaternion
    pub fn exp(&self) -> Quaternion
    {
        // https://en.wikipedia.org/wiki/Quaternion#Functions_of_a_quaternion_variable

        let Quaternion { w, i, j, k } = *self;

        let exp_w = w.exp();
        let axis_norm = (i*i + j*j + k*k).sqrt();

        let (sin,cos) = axis_norm.sin_cos();
        let (sin,cos) = (exp_w * sin / axis_norm, exp_w * cos);

        Quaternion {
            w: cos,
            i: sin * i,
            j: sin * j,
            k: sin * k,
        }

    }

    /// Logarithm of a quaternion
    pub fn log(&self) -> Quaternion
    {
        // https://en.wikipedia.org/wiki/Quaternion#Functions_of_a_quaternion_variable

        let Quaternion { w, i, j, k } = *self;

        let axis_norm = (i*i + j*j + k*k).sqrt();
        let quat_norm = i*i + j*j + k*k + w*w;

        let log  = quat_norm.ln();
        let acos = (w / quat_norm).acos() / axis_norm;

        Quaternion {
            w: log,
            i: acos * i,
            j: acos * j,
            k: acos * k,
        }

    }

}

auto_ops::impl_op_ex!(* |lhs: &Quaternion, rhs: &Quaternion| -> Quaternion {
    Quaternion
    {
        w: lhs.w * rhs.w - lhs.i * rhs.i - lhs.j * rhs.j - lhs.k * rhs.k,
        i: lhs.w * rhs.i + lhs.i * rhs.w + lhs.j * rhs.k - lhs.k * rhs.j,
        j: lhs.w * rhs.j - lhs.i * rhs.k + lhs.j * rhs.w + lhs.k * rhs.i,
        k: lhs.w * rhs.k + lhs.i * rhs.j - lhs.j * rhs.i + lhs.k * rhs.w
    }
});
auto_ops::impl_op_ex_commutative!(* |lhs: &Quaternion, rhs: &Scalar| -> Quaternion {
    Quaternion
    {
        w: lhs.w * rhs,
        i: lhs.i * rhs,
        j: lhs.j * rhs,
        k: lhs.k * rhs
    }
});
auto_ops::impl_op_ex!(*= |lhs: &mut Quaternion, rhs: &Scalar| {
    lhs.w = lhs.w * rhs;
    lhs.i = lhs.i * rhs;
    lhs.j = lhs.j * rhs;
    lhs.k = lhs.k * rhs;
});

auto_ops::impl_op_ex!(/ |lhs: &Quaternion, rhs: &Quaternion| -> Quaternion { lhs * rhs.conj() * (1.0 / rhs.norm().powi(2) ) });
auto_ops::impl_op_ex!(/ |lhs: &Quaternion, rhs: &Scalar| -> Quaternion {
    Quaternion
    {
        w: lhs.w / rhs,
        i: lhs.i / rhs,
        j: lhs.j / rhs,
        k: lhs.k / rhs
    }
});
auto_ops::impl_op_ex!(/ |lhs: &Scalar, rhs: &Quaternion| -> Quaternion { lhs * rhs.conj() * (1.0 / rhs.norm().powi(2) ) });
auto_ops::impl_op_ex!(/= |lhs: &mut Quaternion, rhs: &Scalar| {
    lhs.w /= rhs;
    lhs.i /= rhs;
    lhs.j /= rhs;
    lhs.k /= rhs;
});
