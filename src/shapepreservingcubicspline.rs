use num_traits::{one, zero};

use crate::NewHermite;
use crate::PiecewiseCubicCurve;
use crate::Scalar;
use crate::{fail, Error};

pub trait NewShapePreserving<S: Scalar> {
    fn new_shape_preserving(
        values: impl Into<Vec<S>>,
        grid: impl Into<Vec<S>>,
        closed: bool,
    ) -> Result<PiecewiseCubicCurve<S, S>, Error> {
        let values = values.into();
        let slopes = vec![None; values.len()];
        PiecewiseCubicCurve::new_shape_preserving_with_slopes(values, slopes, grid, closed)
    }
}

impl<S: Scalar> NewShapePreserving<S> for PiecewiseCubicCurve<S, S> {}

pub trait NewShapePreservingWithSlopes<S: Scalar> {
    fn new_shape_preserving_with_slopes(
        values: impl Into<Vec<S>>,
        optional_slopes: impl AsRef<[Option<S>]>,
        grid: impl Into<Vec<S>>,
        closed: bool,
    ) -> Result<PiecewiseCubicCurve<S, S>, Error> {
        let mut values = values.into();
        let optional_slopes = optional_slopes.as_ref();
        let mut grid = grid.into();
        if values.len() < 2 {
            fail!("At least two values are required");
        }
        if values.len() + closed as usize != grid.len() {
            fail!("Number of grid values must be same as values (one more for closed curves)");
        }
        if values.len() != optional_slopes.len() {
            fail!("Number of slopes must be same as values");
        }
        if grid.iter().any(|&x| x.is_nan()) {
            fail!("NaN values are not allowed in grid");
        }
        if closed {
            values.reserve_exact(2);
            grid.reserve_exact(1);
            // Closing the curve:
            values.push(values[0]);
            // One additional (temporary) value to calculate slopes:
            values.push(values[1]);
            grid.push(*grid.last().unwrap() + grid[1] - grid[0]);
        }
        let mut slopes = Vec::<S>::new();
        let two = one::<S>() + one();
        for i in 0..values.len() - 2 {
            let x_1 = values[i];
            let x0 = values[i + 1];
            let x1 = values[i + 2];
            let t_1 = grid[i];
            let t0 = grid[i + 1];
            let t1 = grid[i + 2];
            let left = (x0 - x_1) / (t0 - t_1);
            let right = (x1 - x0) / (t1 - t0);
            let slope = match optional_slopes[(i + 1) % optional_slopes.len()] {
                Some(slope) => {
                    if slope != fix_slope(slope, left, right) {
                        fail!("Slope too steep or wrong sign: {:?}", slope);
                    }
                    slope
                }
                None => fix_slope(
                    ((x0 - x_1) / (t0 - t_1) + (x1 - x0) / (t1 - t0)) / two,
                    left,
                    right,
                ),
            };
            slopes.push(slope); // incoming
            slopes.push(slope); // outgoing
        }
        if closed {
            // Move last (outgoing) tangent to the beginning:
            slopes.rotate_right(1);
            // Remove temporary values:
            values.pop();
            grid.pop();
        } else if optional_slopes.len() == 2 {
            assert!(slopes.is_empty());
            let chord = (values[1] - values[0]) / (grid[1] - grid[0]);
            let one = optional_slopes[0];
            let two = optional_slopes[1];
            slopes.push(check_slope(&one, &two, chord)?);
            slopes.push(check_slope(&two, &one, chord)?);
        } else {
            slopes.insert(
                0,
                check_slope(
                    optional_slopes.first().unwrap(),
                    &slopes.first().copied(),
                    (values[1] - values[0]) / (grid[1] - grid[0]),
                )?,
            );
            slopes.push(check_slope(
                optional_slopes.last().unwrap(),
                &slopes.last().copied(),
                (values[values.len() - 1] - values[values.len() - 2])
                    / (grid[grid.len() - 1] - grid[grid.len() - 2]),
            )?);
        }
        PiecewiseCubicCurve::new_hermite(&values, &slopes, &grid)
    }
}

impl<S: Scalar> NewShapePreservingWithSlopes<S> for PiecewiseCubicCurve<S, S> {}

fn check_slope<S: Scalar>(main: &Option<S>, other: &Option<S>, chord: S) -> Result<S, Error> {
    if let Some(main) = main {
        let main = *main;
        if main != fix_slope(main, chord, chord) {
            fail!("Slope too steep or wrong sign: {:?}", main);
        }
        Ok(main)
    } else if let Some(other) = other {
        Ok(end_slope(*other, chord))
    } else {
        Ok(chord)
    }
}

/// Manipulate the slope to preserve shape.
/// See Dougherty et al. (1989), eq. (4.2).
fn fix_slope<S: Scalar>(slope: S, left: S, right: S) -> S {
    let zero = zero();
    let three = one::<S>() + one() + one();
    if left * right <= zero {
        zero
    } else if right > zero {
        zero.max(slope).min(three * left.abs().min(right.abs()))
    } else {
        zero.min(slope).max(-three * left.abs().min(right.abs()))
    }
}

/// NB: This is a very ad-hoc algorithm meant to minimize the change in slope
/// within the first/last curve segment.  Especially, this should avoid a
/// change from negative to positive acceleration (and vice versa).
/// There might be a better method available!?!
fn end_slope<S: Scalar>(inner_slope: S, chord_slope: S) -> S {
    let zero = zero();
    let two = one::<S>() + one();
    let three = two + one();
    if chord_slope < zero {
        return -end_slope(-inner_slope, -chord_slope);
    }
    assert!(zero <= inner_slope);
    assert!(inner_slope <= three * chord_slope);
    if inner_slope <= chord_slope {
        three * chord_slope - two * inner_slope
    } else {
        (three * chord_slope - inner_slope) / two
    }
}
