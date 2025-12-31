use crate::{complex::Complex, quaternion::Quaternion, dual_quaternion::DualQuaternion, point::Point, direction::Direction};

// Operations between Points and Directions

auto_ops::impl_op_ex!(* |lhs: &Point, rhs: &Point| -> Direction {
    Direction {
        x: lhs.x - rhs.x,
        y: lhs.y - rhs.y,
        z: lhs.z - rhs.z
    }
});

impl From<Point> for Direction
{
    fn from(value: Point) -> Self {
        Direction { x: value.x, y: value.y, z: value.z }
    }
}

impl From<&Point> for Direction
{
    fn from(value: &Point) -> Self {
        Direction { x: value.x, y: value.y, z: value.z }
    }
}

impl From<Direction> for Point
{
    fn from(value: Direction) -> Self {
        Point { x: value.x, y: value.y, z: value.z }
    }
}

impl From<&Direction> for Point
{
    fn from(value: &Direction) -> Self {
        Point { x: value.x, y: value.y, z: value.z }
    }
}

// Operations between Quaternions/Complex numbers

auto_ops::impl_op_ex_commutative!(+ |lhs: &Quaternion, rhs: &Complex| -> Quaternion {
    Quaternion {
        w: lhs.w + rhs.re,
        i: lhs.i + rhs.im,
        j: lhs.j,
        k: lhs.k
    }
});
auto_ops::impl_op_ex!(+= |lhs: &mut Quaternion, rhs: &Complex| {
    lhs.w += rhs.re;
    lhs.i += rhs.im;
});

auto_ops::impl_op_ex!(- |lhs: &Quaternion, rhs: &Complex| -> Quaternion {
    Quaternion {
        w: lhs.w - rhs.re,
        i: lhs.i - rhs.im,
        j: lhs.j,
        k: lhs.k
    }
});
auto_ops::impl_op_ex!(- |lhs: &Complex, rhs: &Quaternion| -> Quaternion {
    Quaternion {
        w: lhs.re - rhs.w,
        i: lhs.im - rhs.i,
        j: -rhs.j,
        k: -rhs.k
    }
});
auto_ops::impl_op_ex!(-= |lhs: &mut Quaternion, rhs: &Complex| {
    lhs.w -= rhs.re;
    lhs.i -= rhs.im;
});

auto_ops::impl_op_ex!(* |lhs: &Quaternion, rhs: &Complex| -> Quaternion {
    Quaternion {
        w:  lhs.w * rhs.re - lhs.i * rhs.im,
        i:  lhs.w * rhs.im + lhs.i * rhs.re,
        j:  lhs.j * rhs.re + lhs.k * rhs.im,
        k: -lhs.j * rhs.im + lhs.k * rhs.re
    }
});
auto_ops::impl_op_ex!(* |lhs: &Complex, rhs: &Quaternion| -> Quaternion {
    Quaternion {
        w: lhs.re * rhs.w - lhs.im * rhs.i,
        i: lhs.re * rhs.i + lhs.im * rhs.w,
        j: lhs.re * rhs.j - lhs.im * rhs.k,
        k: lhs.re * rhs.k + lhs.im * rhs.j
    }
});
auto_ops::impl_op_ex!(*= |lhs: &mut Quaternion, rhs: &Complex| {
    lhs.w =  lhs.w * rhs.re - lhs.i * rhs.im;
    lhs.i =  lhs.w * rhs.im + lhs.i * rhs.re;
    lhs.j =  lhs.j * rhs.re + lhs.k * rhs.im;
    lhs.k = -lhs.j * rhs.im + lhs.k * rhs.re;
});

// These two macros look identical, but they ARE DIFFERENT calculations
auto_ops::impl_op_ex!(/ |lhs: &Quaternion, rhs: &Complex| -> Quaternion { lhs * rhs.conj() * (1.0 / rhs.norm().powi(2) ) });
auto_ops::impl_op_ex!(/ |lhs: &Complex, rhs: &Quaternion| -> Quaternion { lhs * rhs.conj() * (1.0 / rhs.norm().powi(2) ) });
auto_ops::impl_op_ex!(/= |lhs: &mut Quaternion, rhs: &Complex| { *lhs *= rhs.conj() * (1.0 / rhs.norm().powi(2) ); });
