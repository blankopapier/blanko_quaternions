

/// Used instead of f32/f64 to avoid unnecessary generics
#[cfg(not(feature = "use_f64"))]
pub type Scalar = f32;

/// Used instead of f32/f64 to avoid unnecessary generics
#[cfg(feature = "use_f64")]
pub type Scalar = f64;
