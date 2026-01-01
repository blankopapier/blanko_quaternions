use crate::util::Scalar;

#[repr(C)]
#[derive(
    Debug, Clone, Copy, PartialEq, PartialOrd, bytemuck::Pod, bytemuck::Zeroable,
    derive_more::Add, derive_more::AddAssign, derive_more::Sub, derive_more::SubAssign,
    derive_more::Neg, derive_more::From
)]
pub struct Complex
{
    pub re: Scalar,
    pub im: Scalar,
}

impl From<Scalar> for Complex
{
    fn from(value: Scalar) -> Self { Complex { re: value, im: 0.0 } }
}

impl From<&Scalar> for Complex
{
    fn from(value: &Scalar) -> Self { Complex { re: *value, im: 0.0 } }
}

impl std::fmt::Display for Complex
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.re.powi(2) > Scalar::EPSILON {
            write!(f, "{}",  self.re);

            if self.im.powi(2) > Scalar::EPSILON
            {
                write!(f, " + ");
            }
        }

        if self.im.powi(2) > Scalar::EPSILON {
            write!(f, "{}i", self.im);
        }

        write!(f, "")
    }
}

impl Complex
{
    /// Conjugate, i.e. negate the imaginary part
    pub fn conj(&self) -> Self { Self { re: self.re, im: -self.im } }

    /// Norm of this complex number
    pub fn norm(&self) -> Scalar { (self.re*self.re + self.im*self.im).sqrt() }

    /// Normalize this complex number
    pub fn normalized(&self) -> Self { *self * (1.0 / self.norm()) }

    /// May produce invalid numbers if this complex number is 0.0
    pub fn sqrt(&self) -> Self
    {
        //     a+bi = r*exp(ix)
        //  => sqrt(a+bi) = (a+bi)^0.5 = (r*exp(ix))^0.5 = sqrt(r) * exp(ix/2)
        //
        //  r = sqrt(a^2 + b^2);
        //  x = atan2(b,a)

        let x = self.im.atan2(self.re);
        let (sin,cos) = x.sin_cos();

        let r = self.norm();

        Self {
            re: r*cos,
            im: r*sin
        }
    }

    /// Complex exponential function
    pub fn exp(&self) -> Self
    {
        // exp(a+bi) = exp(a)*exp(ib) = exp(a) * (cos(b) + i*sin(b))

        let r = self.re.exp();
        let (sin,cos) = self.im.sin_cos();

        Self {
            re: r*cos,
            im: r*sin
        }
    }

    /// Natural (principal value) logarithm for complex numbers.
    /// This may return invalid numbers if .re <= 0.0
    pub fn log(&self) -> Self
    {
        // https://en.wikipedia.org/wiki/Complex_logarithm

        let im = self.im.atan2(self.re);
        let re = self.norm().ln();

        Self { re, im }
    }

    /// Complex sine function
    pub fn sin(&self) -> Self
    {
        // sin(z) = (exp(iz) - exp(-iz)) / 2i
        //        = (exp(iz) - exp(-iz)) * -0.5i

        // Multiply by i manually
        let x = (Complex { re: -self.im, im: self.re }).exp();
        let y = x - 1.0/x;
        Complex { re: 0.5*y.im, im: -0.5*y.re }
    }

    /// Complex cosine function
    pub fn cos(&self) -> Self
    {
        // cos(z) = (exp(iz) + exp(-iz)) / 2i
        //        = (exp(iz) + exp(-iz)) * -0.5i

        // Multiply by i manually
        let x = (Complex { re: -self.im, im: self.re }).exp();
        let y = x + 1.0/x;
        Complex { re: 0.5*y.im, im: -0.5*y.re }
    }

    /// Complex tangent function
    pub fn tan(&self) -> Self
    {
        // tan(z) = i * (1 - exp(2iz))/(1 + exp(2iz))

        let x = (Complex { re: -2.0 * self.im, im: 2.0 * self.re  }).exp();
        let y = (1.0 - x) / (1.0 + x);
        Complex { re: -y.im, im: y.re }
    }


    /// Raise a complex number to some (real) power.
    /// This may return invalid numbers if this complex number is 0.0
    pub fn powf(&self, f: Scalar) -> Self
    {
        // (a+bi)^n = r^n * exp(i*nx) = r^n * (cos(nx) + i*sin(nx))
        // r = sqrt(a^2 + b^2);
        // x = atan2(b,a);

        let x = self.im.atan2(self.re);
        let (sin,cos) = (f*x).sin_cos();

        let r = self.norm().powf(f);

        Self {
            re: r*cos,
            im: r*sin
        }
    }

    /// Raise a complex number to some integer power.
    /// This may return invalid numbers if this complex number is 0.0
    pub fn powi(&self, i: i32) -> Self
    {
        // (a+bi)^n = r^n * exp(i*nx) = r^n * (cos(nx) + i*sin(nx))
        // r = sqrt(a^2 + b^2);
        // x = atan2(b,a);

        let x = self.im.atan2(self.re);
        let (sin,cos) = ( (i as Scalar) *x).sin_cos();

        let r = self.norm().powi(i);

        Self {
            re: r*cos,
            im: r*sin
        }
    }

    /// Raise a complex number to some complex power.
    /// This may return invalid numbers if this complex number's real part <= 0.0
    pub fn pow(&self, z: Complex) -> Self
    {
        ( z * self.log() ).exp()
    }

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
auto_ops::impl_op_ex_commutative!(* |lhs: &Complex, rhs: &Scalar| -> Complex {
    Complex
    {
        re: lhs.re * rhs,
        im: lhs.im * rhs
    }
});
auto_ops::impl_op_ex!(*= |lhs: &mut Complex, rhs: &Scalar| {
    lhs.re *= rhs;
    lhs.im *= rhs;
});

auto_ops::impl_op_ex!(/ |lhs: &Complex, rhs: &Complex| -> Complex { lhs * rhs.conj() * (1.0 / rhs.norm().powi(2) ) });
auto_ops::impl_op_ex!(/= |lhs: &mut Complex, rhs: &Complex| { *lhs *= rhs.conj() * (1.0 / rhs.norm().powi(2) ) });
auto_ops::impl_op_ex!(/ |lhs: &Complex, rhs: &Scalar| -> Complex {
    Complex
    {
        re: lhs.re / rhs,
        im: lhs.im / rhs
    }
});
auto_ops::impl_op_ex!(/ |lhs: &Scalar, rhs: &Complex| -> Complex { lhs * rhs.conj() * (1.0 / rhs.norm().powi(2) ) });
auto_ops::impl_op_ex!(/= |lhs: &mut Complex, rhs: &Scalar| {
    lhs.re /= rhs;
    lhs.im /= rhs;
});

auto_ops::impl_op_ex_commutative!(+ |lhs: &Complex, rhs: &Scalar| -> Complex {
    Complex
    {
        re: lhs.re + rhs,
        im: lhs.im
    }
});
auto_ops::impl_op_ex!(+= |lhs: &mut Complex, rhs: &Scalar| { lhs.re += rhs });

auto_ops::impl_op_ex!(- |lhs: &Complex, rhs: &Scalar| -> Complex {
    Complex
    {
        re: lhs.re - rhs,
        im: lhs.im
    }
});
auto_ops::impl_op_ex!(- |lhs: &Scalar, rhs: &Complex| -> Complex {
    Complex
    {
        re: lhs - rhs.re,
        im: -rhs.im
    }
});
auto_ops::impl_op_ex!(-= |lhs: &mut Complex, rhs: &Scalar| { lhs.re -= rhs });
