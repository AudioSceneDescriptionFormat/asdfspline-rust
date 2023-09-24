//! Rust implementation of the spline type that's used in the Audio Scene
//! Description Format (ASDF), see
//! <https://AudioSceneDescriptionFormat.readthedocs.io/>.

#![deny(unsafe_code)] // NB: A lot of unsafe code is in the "ffi" sub-crate

use std::ops::{Add, Div, DivAssign, Mul, Sub};

use superslice::Ext; // for slice::upper_bound_by()

pub mod adapters;
pub mod asdfposspline;
pub mod asdfrotspline;
pub mod centripetalkochanekbartelsspline;
pub mod cubichermitespline;
pub mod monotonecubicspline;
pub mod piecewisecubiccurve;
pub mod piecewisemonotonecubicspline;
pub mod quaternion;
pub mod utilities;

pub use crate::asdfposspline::AsdfPosSpline;
pub use crate::asdfrotspline::AsdfRotSpline;
pub use crate::monotonecubicspline::MonotoneCubicSpline;
pub use crate::piecewisecubiccurve::PiecewiseCubicCurve;

use crate::utilities::gauss_legendre13;

/// A trait that is automatically implemented for all types that can be used as positions,
/// polynomial coefficients, tangent vectors etc.
pub trait Vector
where
    Self: Copy,
    Self: Add<Output = Self> + Sub<Output = Self>,
    Self: Mul<f32, Output = Self> + Div<f32, Output = Self>,
    Self: DivAssign<f32>,
{
}

impl<T> Vector for T
where
    Self: Copy,
    Self: Add<Output = Self> + Sub<Output = Self>,
    Self: Mul<f32, Output = Self> + Div<f32, Output = Self>,
    Self: DivAssign<f32>,
{
}

pub trait Spline<Value> {
    fn evaluate(&self, t: f32) -> Value;

    fn grid(&self) -> &[f32];

    /// There must be at least two grid values!
    /// This doesn't work if there are NaNs
    fn clamp_parameter_and_find_index(&self, t: f32) -> (f32, usize) {
        let first = *self.grid().first().unwrap();
        let last = *self.grid().last().unwrap();
        if t < first {
            (first, 0)
        } else if t < last {
            (
                t,
                // NB: This doesn't work if a value is NaN
                self.grid().upper_bound_by(|x| x.partial_cmp(&t).unwrap()) - 1,
            )
        } else {
            (last, self.grid().len() - 2)
        }
    }
}

/// To work around Rust's orphan rules, see <https://blog.mgattozzi.dev/orphan-rules/>
pub trait NormWrapper<U> {
    fn norm(&self) -> f32;
}

pub trait SplineWithVelocity<Value, Velocity>: Spline<Value>
where
    Velocity: Vector,
{
    fn evaluate_velocity(&self, t: f32) -> Velocity;

    fn integrated_speed<U>(&self, index: usize, a: f32, b: f32) -> f32
    where
        Velocity: NormWrapper<U>,
    {
        assert!(a <= b);
        assert!(self.grid()[index] <= a);
        assert!(b <= self.grid()[index + 1]);
        gauss_legendre13(|t| self.evaluate_velocity(t).norm(), a, b)
    }
}
