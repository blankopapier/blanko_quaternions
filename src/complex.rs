// Complex numbers

#[repr(C)]
#[derive(
    Debug, Clone, Copy, PartialEq, PartialOrd, bytemuck::Pod, bytemuck::Zeroable,
    derive_more::Add, derive_more::AddAssign, derive_more::Sub, derive_more::SubAssign,
    derive_more::Neg, derive_more::From
)]
pub struct Complex
{
    pub re: f32,
    pub im: f32,
}

impl From<f32> for Complex
{
    fn from(value: f32) -> Self { Complex { re: value, im: 0.0 } }
}

impl From<&f32> for Complex
{
    fn from(value: &f32) -> Self { Complex { re: *value, im: 0.0 } }
}

impl std::fmt::Display for Complex
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.re.powi(2) > std::f32::EPSILON {
            write!(f, "{}",  self.re);

            if self.im.powi(2) > std::f32::EPSILON
            {
                write!(f, " + ");
            }
        }

        if self.im.powi(2) > std::f32::EPSILON {
            write!(f, "{}i", self.im);
        }

        write!(f, "")
    }
}

impl Complex
{
    pub fn conj(&self) -> Self { Self { re: self.re, im: -self.im } }
    pub fn norm(&self) -> f32 { (self.re*self.re + self.im*self.im).sqrt() }
    pub fn normalized(&self) -> Self { *self * (1.0 / self.norm()) }
}

auto_ops::impl_op_ex!(* |lhs: &Complex, rhs: &Complex| -> Complex {
    Complex
    {
        re: lhs.re * rhs.re - lhs.im * rhs.im,
        im: lhs.im * rhs.re + lhs.re * rhs.im
    }
});
auto_ops::impl_op_ex!(*= |lhs: &mut Complex, rhs: &Complex| {
    lhs.re = lhs.re * rhs.re - lhs.im * rhs.im;
    lhs.im = lhs.im * rhs.re + lhs.re * rhs.im;
});
auto_ops::impl_op_ex_commutative!(* |lhs: &Complex, rhs: &f32| -> Complex {
    Complex
    {
        re: lhs.re * rhs,
        im: lhs.im * rhs
    }
});
auto_ops::impl_op_ex!(*= |lhs: &mut Complex, rhs: &f32| {
    lhs.re *= rhs;
    lhs.im *= rhs;
});

auto_ops::impl_op_ex!(/ |lhs: &Complex, rhs: &Complex| -> Complex { lhs * rhs.conj() * (1.0 / rhs.norm().powi(2) ) });
auto_ops::impl_op_ex!(/= |lhs: &mut Complex, rhs: &Complex| { *lhs *= rhs.conj() * (1.0 / rhs.norm().powi(2) ) });
auto_ops::impl_op_ex!(/ |lhs: &Complex, rhs: &f32| -> Complex {
    Complex
    {
        re: lhs.re / rhs,
        im: lhs.im / rhs
    }
});
auto_ops::impl_op_ex!(/ |lhs: &f32, rhs: &Complex| -> Complex { lhs * rhs.conj() * (1.0 / rhs.norm().powi(2) ) });
auto_ops::impl_op_ex!(/= |lhs: &mut Complex, rhs: &f32| {
    lhs.re /= rhs;
    lhs.im /= rhs;
});

auto_ops::impl_op_ex_commutative!(+ |lhs: &Complex, rhs: &f32| -> Complex {
    Complex
    {
        re: lhs.re + rhs,
        im: lhs.im
    }
});
auto_ops::impl_op_ex!(+= |lhs: &mut Complex, rhs: &f32| { lhs.re += rhs });

auto_ops::impl_op_ex!(- |lhs: &Complex, rhs: &f32| -> Complex {
    Complex
    {
        re: lhs.re - rhs,
        im: lhs.im
    }
});
auto_ops::impl_op_ex!(- |lhs: &f32, rhs: &Complex| -> Complex {
    Complex
    {
        re: lhs - rhs.re,
        im: -rhs.im
    }
});
auto_ops::impl_op_ex!(-= |lhs: &mut Complex, rhs: &f32| { lhs.re -= rhs });
