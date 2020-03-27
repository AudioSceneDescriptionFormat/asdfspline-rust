/*!
Rust implementation of the spline type that's used in the Audio Scene
Description Format (ASDF), see
<https://AudioSceneDescriptionFormat.readthedocs.io/>.

# Requirements

* Rust compiler, Cargo (<https://rustup.rs/>)

The required Rust packages (a.k.a. "crates") are listed in the file
`Cargo.toml`.

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

mod asdfspline;
mod centripetalkochanekbartelsspline;
mod cubichermitespline;
mod error;
mod monotonecubicspline;
mod piecewisecubiccurve;
mod shapepreservingcubicspline;
pub mod utilities;

pub use crate::asdfspline::AsdfSpline;
pub use crate::centripetalkochanekbartelsspline::NewCentripetalKochanekBartels;
pub use crate::cubichermitespline::NewHermite;
pub use crate::error::Error;
pub use crate::monotonecubicspline::MonotoneCubicSpline;
pub use crate::piecewisecubiccurve::PiecewiseCubicCurve;
pub use shapepreservingcubicspline::{NewShapePreserving, NewShapePreservingWithSlopes};

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
