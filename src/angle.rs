//! Some people prefer giving angles in radians while others prefer degrees.
//! This type solves this problem by containing both radians and degrees, while
//! letting the user decide which preference they have.
//!
//! `Angle::new(angle: f32)` will take in radians by default.
//! You can change this to use degrees by default by enabling the `"angle_new_degrees"` feature

use crate::util::Scalar;

#[cfg(not(feature = "use_f64"))]
use std::f32::consts::{PI, TAU};

#[cfg(feature = "use_f64")]
use std::f64::consts::{PI, TAU};

/// Positive angles are counter-clockwise (ccw)

#[repr(C)]
#[derive(
    Debug, Clone, Copy, PartialEq, PartialOrd, bytemuck::Pod, bytemuck::Zeroable,
    derive_more::Add, derive_more::AddAssign, derive_more::Sub, derive_more::SubAssign,
    derive_more::Rem, derive_more::RemAssign, derive_more::Neg
)]
pub struct Angle
{
    rad: Scalar,
    deg: Scalar,
}

impl std::fmt::Display for Angle
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Angle {{ {}° / {}π }}", self.deg, self.rad)
    }
}

auto_ops::impl_op_ex_commutative!(* |lhs: &Angle, rhs: &Scalar| -> Angle {
    Angle { rad: lhs.rad * rhs, deg: lhs.deg* rhs }
});
auto_ops::impl_op_ex!(*= |lhs: &mut Angle, rhs: &Scalar| {
    lhs.rad *= rhs;
    lhs.deg *= rhs;
});

auto_ops::impl_op_ex!(/ |lhs: &Angle, rhs: &Scalar| -> Angle {
    Angle { rad: lhs.rad / rhs, deg: lhs.deg / rhs }
});
auto_ops::impl_op_ex!(/= |lhs: &mut Angle, rhs: &Scalar| {
    lhs.rad /= rhs;
    lhs.deg /= rhs;
});

impl Angle
{
    /// 360°
    pub const FULL   : Angle = Angle { rad: TAU,    deg: 360.0 };

    /// 180°
    pub const HALF   : Angle = Angle { rad: PI,     deg: 180.0 };

    /// 90°
    pub const QUARTER: Angle = Angle { rad: PI/2.0, deg:  90.0 };

    /// 45°
    pub const EIGTH:   Angle = Angle { rad: PI/4.0, deg:  45.0 };

    /// 0°
    pub const ZERO   : Angle = Angle { rad: 0.0,    deg:   0.0 };

    /// Create Angle from degrees
    #[cfg(feature = "angle_new_degrees")]
    pub fn new(angle: Scalar) -> Self
    {
        Self::degrees(angle)
    }

    /// Create Angle from radians
    #[cfg(not(feature = "angle_new_degrees"))]
    pub fn new(angle: Scalar) -> Self
    {
        Self::radians(angle)
    }

    /// Create Angle from radians
    pub fn radians(angle: Scalar) -> Self
    {
        Self { rad: angle, deg: angle*(180.0/PI) }
    }

    /// Create Angle from degrees
    pub fn degrees(angle: Scalar) -> Self
    {
        Self { deg: angle, rad: angle*(PI/180.0) }
    }

    /// Get this angle in radians
    pub fn rad(&self) -> Scalar
    {
        self.rad
    }

    /// Get this angle in degrees
    pub fn deg(&self) -> Scalar
    {
        self.deg
    }

    /// Minimum of two angles
    pub fn min(&self, angle: Angle) -> Angle
    {
        if self.deg < angle.deg { *self } else { angle }
    }

    /// Maximum of two angles
    pub fn max(&self, angle: Angle) -> Angle
    {
        if self.deg > angle.deg { *self } else { angle }
    }

    /// Clamp this Angle between two angles
    pub fn clamp(&self, min: Angle, max: Angle) -> Angle
    {
        if self.deg > max.deg { max }
        else if self.deg < min.deg { min }
        else { *self }
    }

    /// Ignore sign of angle.
    pub fn abs(&self) -> Self
    {
        Self { rad: self.rad.abs(), deg: self.deg.abs() }
    }

    pub fn signum(&self) -> Scalar
    {
        self.deg.signum()
    }

    /// Round to the nearest smaller multiple of some Angle
    pub fn floor(&self, angle: Angle) -> Self
    {
        let n = (self.deg / angle.deg).trunc();

        angle * n
    }

    /// Round to the nearest greater multiple of some Angle
    pub fn ceil(&self, angle: Angle) -> Self
    {
        let a = self.deg / angle.deg;
        let n = if a < 0.0 {
            a.floor()
        } else {
            a.ceil()
        };

        angle * n
    }

    /// Round to the nearest multiple of some Angle
    pub fn round(&self, angle: Angle) -> Self
    {
        let n = (self.deg / angle.deg).round();
        angle * n
    }

    /// Will this Angle to its positive value in [0°,360°).
    /// If this Angle would be -20°, then this method will return 340°
    pub fn corrected(&self) -> Self
    {
        let modulo = self.deg % 360.0;
        let angle  = if self.deg < 0.0 { 360.0 + modulo } else { modulo };
        Self::radians(angle)
    }

    /// Will convert this Angle to its positive value without clamping to [0°,360°).
    /// If this Angle would be -380°, then this method will return 700°
    pub fn sign_corrected(&self) -> Self
    {
        let angle = if self.deg < 0.0
        {
            let n_deg = (self.deg / 360.0).floor() * -360.0;
            n_deg + (self.deg % 360.0)
        }
        else
        {
            self.deg
        };

        Self::degrees(angle)
    }

    /// Will return the correct angle in (-360°,360°).
    /// If this Angle would be -20°, then this method will return 340°
    pub fn range_corrected(&self) -> Self
    {
        Self::degrees(self.deg % 360.0)
    }

    pub fn sin(&self) -> Scalar { self.rad.sin() }
    pub fn cos(&self) -> Scalar { self.rad.cos() }
    pub fn tan(&self) -> Scalar { self.rad.tan() }

    pub fn sin_cos(&self) -> (Scalar,Scalar) { self.rad.sin_cos() }
}
