//! Dual numbers are like complex numbers, but instead of i^2 = -1 we define E^2 = 0.
//! They can be used for some niche applications and mechanical stuff

use crate::util::Scalar;

#[repr(C)]
#[derive(
    Debug, Clone, Copy, PartialEq, PartialOrd, bytemuck::Pod, bytemuck::Zeroable,
    derive_more::Add, derive_more::AddAssign, derive_more::Sub, derive_more::SubAssign,
    derive_more::Neg, derive_more::From
)]
pub struct DualNumber
{
    pub re: Scalar,
    pub du: Scalar,
}

impl From<Scalar> for DualNumber
{
    fn from(value: Scalar) -> Self { DualNumber { re: value, du: 0.0 } }
}

impl From<&Scalar> for DualNumber
{
    fn from(value: &Scalar) -> Self { DualNumber { re: *value, du: 0.0 } }
}

impl std::fmt::Display for DualNumber
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.re.powi(2) > Scalar::EPSILON {
            write!(f, "{}",  self.re);

            if self.du.powi(2) > Scalar::EPSILON
            {
                write!(f, " + ");
            }
        }

        if self.du.powi(2) > Scalar::EPSILON {
            write!(f, "{}i", self.du);
        }

        write!(f, "")
    }
}

impl DualNumber
{
    /// Conjugate, i.e. negate the dual part
    pub fn conj(&self) -> Self { Self { re: self.re, du: -self.du } }

    /// This is the "natural" norm defined via conjugation (a*a.conj).
    /// If this returns 0.0, you can not know whether or not the DualNumber is actually zero (dual part may be non-zero)
    pub fn seminorm(&self) -> Scalar { (self.re*self.re).sqrt() }

    /// Normalizes this DualNumber by the seminorm
    pub fn seminormalized(&self) -> Self { *self * (1.0 / self.seminorm()) }

    /// This is an "artificial" norm equal to the Euclidean Norm.
    /// The DualNumber will be 0, if this norm returns 0.0
    pub fn norm(&self) -> Scalar { (self.re*self.re + self.du*self.du).sqrt() }

    /// Normalizes this DualNumber by the Euclidean Norm
    pub fn normalized(&self) -> Self { *self * (1.0 / self.norm()) }

    /// This may return invalid numbers if .re <= 0.0
    pub fn sqrt(&self) -> Self
    {
        //    sqrt(a+bE) = (A+BE)
        // => a+bE = (A+BE)^2 = A^2 + 2ABE
        // => a = A^2; b = 2AB
        // => A = sqrt(a) => B = b/2A = b/2sqrt(a)

        let s = self.re.sqrt();
        Self { re: s, du: self.du/(2.0*s) }
    }

    pub fn exp(&self) -> Self
    {
        // (a+bE)^n = ... (starting with n=0)
        // 1
        // (a+bE)
        // (a²+2abE)
        // (a³+3a²bE)
        // (a^4 + 4a³bE)
        // ...
        // (a^n + n*a^(n-1)*bE)
        //
        // => exp(a+bE) = exp(a) + bE * exp(a) = exp(a)*(1+bE)

        let exp = self.re.exp();
        Self {
            re: exp,
            du: self.du * exp
        }
    }

    /// Natural logarithm for DualNumbers.
    /// This may return invalid numbers if .re <= 0.0
    pub fn log(&self) -> Self
    {
        //    log(a+bE) = A+BE
        //  => exp(log(a+bE)) = exp(A+BE)
        // <=> a+bE = exp(A) + E*exp(A)*B
        //  => a = exp(A); b = exp(A)*B
        //  => A = log(a); B = b/a
        //  => log(a+bE) = log(a) + E*b/a
        Self {
            re: self.re.ln(),
            du: self.du/self.re
        }
    }

    /// Dual number sine function
    pub fn sin(&self) -> Self
    {
        // z = a + bE
        // sin(z) = Σ (-1)^k * z^(2k+1) / (2k+1)!    k = 0,1,2,3...
        //        = sin(a) + bE*cos(a)

        DualNumber { re: self.re.sin(), du: self.re.cos() * self.du }
    }

    /// Dual number cosine function
    pub fn cos(&self) -> Self
    {
        // For (real) analytic functions, you can extend them via Taylor series.
        // let z = a + bE.
        // Evaluate f(a).
        //
        // Taylor Series:
        // f(a+bE) = Σ f^(k)(a) * (a+bE - a)^k / k!    k = 0,1,2,3...
        //         = Σ f^(k)(a) * (bE)^k / k!          E^k = 0, when k >= 1
        //         = f(a) + f'(a)bE

        DualNumber { re: self.re.cos(), du: -self.re.sin() * self.du }
    }

    /// Dual number tangent function.
    /// May produce invalid numbers when .re = (2k+1)*PI/2 (odd multiple of Pi/2)
    pub fn tan(&self) -> Self
    {
        // f(a+bE) = f(a) + f'(a)bE
        // tan' = 1 / cos^2
        DualNumber { re: self.re.tan(), du: self.du / self.re.cos().powi(2) }
    }

    /// Raise a DualNumber to some (real) power.
    /// This may return invalid numbers if .re <= 0.0
    pub fn powf(&self, f: Scalar) -> Self
    {
        ( f * self.log() ).exp()
    }

    /// Raise a DualNumber to some integer power.
    /// This may return invalid numbers if .re <= 0.0
    pub fn powi(&self, i: i32) -> Self
    {
        // (a+bE)^n = ... (starting with n=0)
        // 1
        // (a+bE)
        // (a²+2abE)
        // (a³+3a²bE)
        // (a^4 + 4a³bE)
        // ...
        // (a^n + n*a^(n-1)*bE)

        // for n < 0, just divide 1 by (a+bE)^n

        let p = self.re.powi( (i-1).max(0) );
        let d = DualNumber {
            re: p*self.re,
            du: (i as Scalar)*p*self.du
        };

        if i < 0 { 1.0 / d } else { d }
    }
}



auto_ops::impl_op_ex!(* |lhs: &DualNumber, rhs: &DualNumber| -> DualNumber {
    DualNumber {
        re: lhs.re * rhs.re,
        du: lhs.re * rhs.du + lhs.du * rhs.re,
    }
});
auto_ops::impl_op_ex!(*= |lhs: &mut DualNumber, rhs: &DualNumber| {
    lhs.re = lhs.re * rhs.re;
    lhs.du = lhs.du * rhs.re + lhs.re * rhs.du;
});
auto_ops::impl_op_ex_commutative!(* |lhs: &DualNumber, rhs: &Scalar| -> DualNumber {
    DualNumber
    {
        re: lhs.re * rhs,
        du: lhs.du * rhs
    }
});
auto_ops::impl_op_ex!(*= |lhs: &mut DualNumber, rhs: &Scalar| {
    lhs.re *= rhs;
    lhs.du *= rhs;
});

auto_ops::impl_op_ex!(/ |lhs: &DualNumber, rhs: &DualNumber| -> DualNumber { lhs * rhs.conj() * (1.0 / rhs.seminorm().powi(2) ) });
auto_ops::impl_op_ex!(/= |lhs: &mut DualNumber, rhs: &DualNumber| { *lhs *= rhs.conj() * (1.0 / rhs.seminorm().powi(2) ) });
auto_ops::impl_op_ex!(/ |lhs: &DualNumber, rhs: &Scalar| -> DualNumber {
    DualNumber
    {
        re: lhs.re / rhs,
        du: lhs.du / rhs
    }
});
auto_ops::impl_op_ex!(/ |lhs: &Scalar, rhs: &DualNumber| -> DualNumber { lhs * rhs.conj() * (1.0 / rhs.seminorm().powi(2) ) });
auto_ops::impl_op_ex!(/= |lhs: &mut DualNumber, rhs: &Scalar| {
    lhs.re /= rhs;
    lhs.du /= rhs;
});

auto_ops::impl_op_ex_commutative!(+ |lhs: &DualNumber, rhs: &Scalar| -> DualNumber {
    DualNumber
    {
        re: lhs.re + rhs,
        du: lhs.du
    }
});
auto_ops::impl_op_ex!(+= |lhs: &mut DualNumber, rhs: &Scalar| { lhs.re += rhs });

auto_ops::impl_op_ex!(- |lhs: &DualNumber, rhs: &Scalar| -> DualNumber {
    DualNumber
    {
        re: lhs.re - rhs,
        du: lhs.du
    }
});
auto_ops::impl_op_ex!(- |lhs: &Scalar, rhs: &DualNumber| -> DualNumber {
    DualNumber
    {
        re: lhs - rhs.re,
        du: -rhs.du
    }
});
auto_ops::impl_op_ex!(-= |lhs: &mut DualNumber, rhs: &Scalar| { lhs.re -= rhs });
