#[repr(C)]
#[derive(
    Debug, Clone, Copy, PartialEq, bytemuck::Pod, bytemuck::Zeroable,
    derive_more::Add, derive_more::AddAssign, derive_more::Sub, derive_more::SubAssign,
    derive_more::Neg, derive_more::From
)]
pub struct Direction
{
    pub x: f32,
    pub y: f32,
    pub z: f32
}

impl Direction
{
    pub fn norm(&self) -> f32
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

auto_ops::impl_op_ex_commutative!(* |lhs: &Direction, rhs: &f32| -> Direction {
    Direction {
        x: lhs.x * rhs,
        y: lhs.y * rhs,
        z: lhs.z * rhs
    }
});

#[macro_export]
macro_rules! direction {
    ($x: expr, $y: expr, $z: expr) => {
        Direction { x: ($x) as f32, y: ($y) as f32, z: ($z) as f32 }
    };
}
