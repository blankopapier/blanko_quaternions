
#[repr(C)]
#[derive(
    Debug, Clone, Copy, PartialEq, bytemuck::Pod, bytemuck::Zeroable,
    derive_more::Neg, derive_more::From
)]
pub struct Point
{
    pub x: f32,
    pub y: f32,
    pub z: f32
}

impl Point
{
    pub fn norm(&self) -> f32
    {
        (self.x*self.x + self.y*self.y + self.z*self.z).sqrt()
    }

    pub fn normalize(&self) -> Self
    {
        self * (1.0/self.norm())
    }
}

auto_ops::impl_op_ex_commutative!(* |lhs: &Point, rhs: &f32| -> Point {
    Point {
        x: lhs.x * rhs,
        y: lhs.y * rhs,
        z: lhs.z * rhs
    }
});

#[macro_export]
macro_rules! point {
    ($x: expr, $y: expr, $z: expr) => {
        Point { x: ($x) as f32, y: ($y) as f32, z: ($z) as f32 }
    };
}
