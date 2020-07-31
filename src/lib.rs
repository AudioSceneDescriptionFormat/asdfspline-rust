/*!
Rust implementation of the spline type that's used in the Audio Scene
Description Format (ASDF), see
<https://AudioSceneDescriptionFormat.readthedocs.io/>.

# Requirements

* Rust compiler, Cargo (<https://rustup.rs/>)

The required Rust packages (a.k.a. "crates") are listed in the file
`Cargo.toml`.

# Tests

```text
cargo test --all
```

There are further tests (using Python) in the `python/` directory.

# API Documentation

Run `cargo doc --all` in the main directory to create the documentation.
The generated HTML documentation can be accessed via
[target/doc/asdfspline/index.html](index.html) and
[target/doc/asdfspline_ffi/index.html](../asdfspline_ffi/index.html).

# Updating `README.md`

Using [cargo-readme](https://github.com/livioribeiro/cargo-readme) (`cargo install cargo-readme`):

```text
cargo readme -o README.md
```
*/
use std::ops::{Add, Div, DivAssign, Mul, Sub};

use superslice::Ext; // for slice::upper_bound_by()

pub mod adapters;
pub mod asdfposspline;
pub mod centripetalkochanekbartelsspline;
pub mod cubichermitespline;
pub mod monotonecubicspline;
pub mod piecewisecubiccurve;
pub mod quaternion;
pub mod shapepreservingcubicspline;
pub mod utilities;

pub use crate::asdfposspline::AsdfPosSpline;
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

/// To work around Rust's orphan rules, see https://blog.mgattozzi.dev/orphan-rules/
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
