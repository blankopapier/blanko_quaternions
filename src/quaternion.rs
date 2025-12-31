// Quaternions

pub use crate::angle::Angle;
pub use crate::point::Point;
pub use crate::direction::Direction;

#[repr(C)]
#[derive(
    Debug, Clone, Copy, PartialEq, bytemuck::Pod, bytemuck::Zeroable,
    derive_more::Add, derive_more::AddAssign, derive_more::Sub, derive_more::SubAssign,
    derive_more::Neg, derive_more::From
)]
pub struct Quaternion
{
    pub w: f32,
    pub i: f32,
    pub j: f32,
    pub k: f32,
}

impl std::fmt::Display for Quaternion
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.w.powi(2) > std::f32::EPSILON {
            write!(f, "{}",  self.w);

            if self.i.powi(2) > std::f32::EPSILON ||
                self.j.powi(2) > std::f32::EPSILON ||
                self.k.powi(2) > std::f32::EPSILON
            {
                write!(f, " + ");
            }
        }

        if self.i.powi(2) > std::f32::EPSILON {
            write!(f, "{}i", self.i);

            if self.j.powi(2) > std::f32::EPSILON ||
                self.k.powi(2) > std::f32::EPSILON
            {
                write!(f, " + ");
            }
        }

        if self.j.powi(2) > std::f32::EPSILON {
            write!(f, "{}j", self.j);

            if self.k.powi(2) > std::f32::EPSILON
            {
                write!(f, " + ");
            }
        }

        if self.k.powi(2) > std::f32::EPSILON {
            write!(f, "{}k", self.k);
        }

        write!(f, "")
    }
}

impl Quaternion
{
    pub fn conj(&self) -> Self { Self { w: self.w, i: -self.i, j: -self.j, k: -self.k } }
    pub fn norm(&self) -> f32 { (self.w*self.w + self.i*self.i + self.j*self.j + self.k*self.k).sqrt() }
    pub fn normalized(&self) -> Self { *self * (1.0 / self.norm()) }

    pub fn from_angle_axis(angle: Angle, x: f32, y: f32, z: f32) -> Self
    {
        let mut q = Self { w: 0.0, i: x, j: y, k: z }.normalized();
        let (sin,cos) = (angle*0.5).sin_cos();

        q *= sin;
        q.w = cos;

        q
    }

    pub fn transform_point(&self, point: &Point) -> Point
    {
        self.transform_direction(&point.into()).into()
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

auto_ops::impl_op_ex!(* |lhs: &Quaternion, rhs: &Quaternion| -> Quaternion {
    Quaternion
    {
        w: lhs.w * rhs.w - lhs.i * rhs.i - lhs.j * rhs.j - lhs.k * rhs.k,
        i: lhs.w * rhs.i + lhs.i * rhs.w + lhs.j * rhs.k - lhs.k * rhs.j,
        j: lhs.w * rhs.j - lhs.i * rhs.k + lhs.j * rhs.w + lhs.k * rhs.i,
        k: lhs.w * rhs.k + lhs.i * rhs.j - lhs.j * rhs.i + lhs.k * rhs.w
    }
});
auto_ops::impl_op_ex_commutative!(* |lhs: &Quaternion, rhs: &f32| -> Quaternion {
    Quaternion
    {
        w: lhs.w * rhs,
        i: lhs.i * rhs,
        j: lhs.j * rhs,
        k: lhs.k * rhs
    }
});
auto_ops::impl_op_ex!(*= |lhs: &mut Quaternion, rhs: &f32| {
    lhs.w = lhs.w * rhs;
    lhs.i = lhs.i * rhs;
    lhs.j = lhs.j * rhs;
    lhs.k = lhs.k * rhs;
});

auto_ops::impl_op_ex!(/ |lhs: &Quaternion, rhs: &Quaternion| -> Quaternion { lhs * rhs.conj() * (1.0 / rhs.norm().powi(2) ) });
auto_ops::impl_op_ex!(/ |lhs: &Quaternion, rhs: &f32| -> Quaternion {
    Quaternion
    {
        w: lhs.w / rhs,
        i: lhs.i / rhs,
        j: lhs.j / rhs,
        k: lhs.k / rhs
    }
});
auto_ops::impl_op_ex!(/ |lhs: &f32, rhs: &Quaternion| -> Quaternion { lhs * rhs.conj() * (1.0 / rhs.norm().powi(2) ) });
auto_ops::impl_op_ex!(/= |lhs: &mut Quaternion, rhs: &f32| {
    lhs.w /= rhs;
    lhs.i /= rhs;
    lhs.j /= rhs;
    lhs.k /= rhs;
});
