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
use std::fmt::Debug;
use std::ops::{Add, Div, DivAssign, Mul, Sub};

use num_traits::{Float, FromPrimitive, NumAssign};

pub mod adapters;
pub mod asdfspline;
pub mod centripetalkochanekbartelsspline;
pub mod cubichermitespline;
pub mod monotonecubicspline;
pub mod piecewisecubiccurve;
pub mod shapepreservingcubicspline;
pub mod utilities;

pub use crate::asdfspline::AsdfPosSpline;
pub use crate::monotonecubicspline::MonotoneCubicSpline;
pub use crate::piecewisecubiccurve::PiecewiseCubicCurve;

use crate::utilities::gauss_legendre13;

/// A trait that is automatically implemented for all types that can be used as scalars,
/// e.g. time values.
pub trait Scalar: Float + NumAssign + FromPrimitive + Debug {}

impl<T: Float + NumAssign + FromPrimitive + Debug> Scalar for T {}

/// A trait that is automatically implemented for all types that can be used as positions,
/// polynomial coefficients, tangent vectors etc.
pub trait Vector<S>
where
    S: Scalar,
    Self: Copy,
    Self: Add<Output = Self> + Sub<Output = Self>,
    Self: Mul<S, Output = Self> + Div<S, Output = Self>,
    Self: DivAssign<S>,
{
}

impl<S, T> Vector<S> for T
where
    S: Scalar,
    Self: Copy,
    Self: Add<Output = Self> + Sub<Output = Self>,
    Self: Mul<S, Output = Self> + Div<S, Output = Self>,
    Self: DivAssign<S>,
{
}

pub trait Spline<S, Output>
where
    S: Scalar,
{
    fn evaluate(&self, t: S) -> Output;
    fn grid(&self) -> &[S];
}

/// To work around Rust's orphan rules, see https://blog.mgattozzi.dev/orphan-rules/
pub trait NormWrapper<Dummy> {
    type Norm;
    fn norm(&self) -> Self::Norm;
}

pub trait SplineWithVelocity<S, Output, Velocity>: Spline<S, Output>
where
    S: Scalar,
    Velocity: Vector<S>,
{
    fn evaluate_velocity(&self, t: S) -> Velocity;

    fn integrated_speed<Dummy>(&self, index: usize, a: S, b: S) -> S
    where
        Velocity: NormWrapper<Dummy, Norm = S>,
    {
        assert!(a <= b);
        assert!(self.grid()[index] <= a);
        assert!(b <= self.grid()[index + 1]);
        gauss_legendre13(|t| self.evaluate_velocity(t).norm(), a, b)
    }
}
