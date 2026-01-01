//! This crate contains types and functionality related to:
//!
//! * Complex numbers (`Complex`)
//! * Dual numbers (`DualNumber`)
//! * Quaternions (`Quaternion`)
//! * Dual quaternions (`DualQuaternion`)
//! * Angles (`Angle`)
//!
//! <div class="warning">
//! This crate is still in development, but usable.
//! Therefore it may still change quite a bit
//! </div>
//!
//! # Cargo features
//! * `angle_new_degrees` will make `Angle::new(angle)` use degrees as input (disabled by default)
//! * `use_f64` will use f64 as scalar type for components instead of f32 (disabled by default)



pub mod angle;

pub mod complex;
pub mod dual_numbers;

pub mod quaternion;
pub mod dual_quaternion;

mod util;
mod vector3;
