

#[repr(C)]
#[derive(
    Debug, Clone, Copy, PartialEq, PartialOrd, bytemuck::Pod, bytemuck::Zeroable,
    derive_more::Add, derive_more::AddAssign, derive_more::Sub, derive_more::SubAssign,
    derive_more::Neg, derive_more::From
)]
pub struct DualNumber
{
    pub re: f32,
    pub du: f32,
}

impl From<f32> for DualNumber
{
    fn from(value: f32) -> Self { DualNumber { re: value, du: 0.0 } }
}

impl From<&f32> for DualNumber
{
    fn from(value: &f32) -> Self { DualNumber { re: *value, du: 0.0 } }
}

impl std::fmt::Display for DualNumber
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.re.powi(2) > std::f32::EPSILON {
            write!(f, "{}",  self.re);

            if self.du.powi(2) > std::f32::EPSILON
            {
                write!(f, " + ");
            }
        }

        if self.du.powi(2) > std::f32::EPSILON {
            write!(f, "{}i", self.du);
        }

        write!(f, "")
    }
}

impl DualNumber
{
    pub fn conj(&self) -> Self { Self { re: self.re, du: -self.du } }

    /// This is the "natural" norm defined via conjugation (a*a.conj).
    /// If this returns 0.0, you can not know whether or not the DualNumber is actually zero (dual part may be non-zero)
    pub fn seminorm(&self) -> f32 { (self.re*self.re).sqrt() }

    /// Normalizes this DualNumber by the seminorm
    pub fn seminormalized(&self) -> Self { *self * (1.0 / self.seminorm()) }

    /// This is an "artificial" norm equal to the Euclidean Norm.
    /// The DualNumber will be 0, if this norm returns 0.0
    pub fn norm(&self) -> f32 { (self.re*self.re + self.du*self.du).sqrt() }

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

    /// Raise a DualNumber to some (real) power.
    /// This may return invalid numbers if .re <= 0.0
    pub fn powf(&self, f: f32) -> Self
    {
        ( f * self.log() ).exp()
    }

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
            du: (i as f32)*p*self.du
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
auto_ops::impl_op_ex_commutative!(* |lhs: &DualNumber, rhs: &f32| -> DualNumber {
    DualNumber
    {
        re: lhs.re * rhs,
        du: lhs.du * rhs
    }
});
auto_ops::impl_op_ex!(*= |lhs: &mut DualNumber, rhs: &f32| {
    lhs.re *= rhs;
    lhs.du *= rhs;
});

auto_ops::impl_op_ex!(/ |lhs: &DualNumber, rhs: &DualNumber| -> DualNumber { lhs * rhs.conj() * (1.0 / rhs.seminorm().powi(2) ) });
auto_ops::impl_op_ex!(/= |lhs: &mut DualNumber, rhs: &DualNumber| { *lhs *= rhs.conj() * (1.0 / rhs.seminorm().powi(2) ) });
auto_ops::impl_op_ex!(/ |lhs: &DualNumber, rhs: &f32| -> DualNumber {
    DualNumber
    {
        re: lhs.re / rhs,
        du: lhs.du / rhs
    }
});
auto_ops::impl_op_ex!(/ |lhs: &f32, rhs: &DualNumber| -> DualNumber { lhs * rhs.conj() * (1.0 / rhs.seminorm().powi(2) ) });
auto_ops::impl_op_ex!(/= |lhs: &mut DualNumber, rhs: &f32| {
    lhs.re /= rhs;
    lhs.du /= rhs;
});

auto_ops::impl_op_ex_commutative!(+ |lhs: &DualNumber, rhs: &f32| -> DualNumber {
    DualNumber
    {
        re: lhs.re + rhs,
        du: lhs.du
    }
});
auto_ops::impl_op_ex!(+= |lhs: &mut DualNumber, rhs: &f32| { lhs.re += rhs });

auto_ops::impl_op_ex!(- |lhs: &DualNumber, rhs: &f32| -> DualNumber {
    DualNumber
    {
        re: lhs.re - rhs,
        du: lhs.du
    }
});
auto_ops::impl_op_ex!(- |lhs: &f32, rhs: &DualNumber| -> DualNumber {
    DualNumber
    {
        re: lhs - rhs.re,
        du: -rhs.du
    }
});
auto_ops::impl_op_ex!(-= |lhs: &mut DualNumber, rhs: &f32| { lhs.re -= rhs });
