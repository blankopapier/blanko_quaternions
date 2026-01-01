use crate::util::Scalar;

/// Used internally for easier manipulation of vectors.
/// Internally only to keep this library simple and compatible with other, probably better, lin-alg crates
#[repr(C)]
#[derive(
    Debug, Clone, Copy, PartialEq, bytemuck::Pod, bytemuck::Zeroable,
    derive_more::Add, derive_more::AddAssign, derive_more::Sub, derive_more::SubAssign,
    derive_more::Neg, derive_more::From
)]
pub struct Vector3
{
    pub x: Scalar,
    pub y: Scalar,
    pub z: Scalar
}

impl Vector3
{
    pub fn norm(&self) -> Scalar
    {
        ( self.x*self.x + self.y*self.y + self.z*self.z).sqrt()
    }

    pub fn normalize(&self) -> Self
    {
        self * (1.0/self.norm())
    }

    pub fn cross(&self, other: &Self) -> Self
    {
        Self {
            x: self.y * other.z - other.y * self.z,
            y: self.z * other.x - other.z * self.x,
            z: self.x * other.y - other.x * self.y
        }
    }
}

auto_ops::impl_op_ex_commutative!(* |lhs: &Vector3, rhs: &Scalar| -> Vector3 {
    Vector3 {
        x: lhs.x * rhs,
        y: lhs.y * rhs,
        z: lhs.z * rhs
    }
});

impl From<Vector3> for [Scalar;3]
{
    fn from(value: Vector3) -> Self {
        [value.x, value.y, value.z]
    }
}
