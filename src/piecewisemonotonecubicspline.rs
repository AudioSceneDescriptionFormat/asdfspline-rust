use crate::utilities::{check_grid, GridError};
use crate::PiecewiseCubicCurve;

#[derive(thiserror::Error, Debug)]
pub enum PiecewiseMonotoneError {
    #[error("there must be at least two values")]
    LessThanTwoValues,
    #[error("length of grid ({grid}) must be {} number of values ({values})", if *.closed {
        "one more than"
    } else {
        "the same as"
    })]
    GridVsValues {
        grid: usize,
        values: usize,
        closed: bool,
    },
    #[error(transparent)]
    FromGridError(#[from] GridError),
}

#[derive(thiserror::Error, Debug)]
pub enum PiecewiseMonotoneWithSlopesError {
    #[error(transparent)]
    FromPiecewiseMonotoneError(#[from] PiecewiseMonotoneError),
    #[error("number of slopes ({slopes}) must be same as number of values ({values})")]
    SlopesVsValues { slopes: usize, values: usize },
    #[error("slope at index {index} too steep ({slope:?}; maximum: {maximum:?})")]
    SlopeTooSteep {
        index: usize,
        slope: f32,
        maximum: f32,
    },
    #[error("slope at index {index} has wrong sign ({slope:?})")]
    SlopeWrongSign { index: usize, slope: f32 },
}

impl PiecewiseCubicCurve<f32> {
    pub fn new_piecewise_monotone(
        values: impl Into<Vec<f32>>,
        grid: impl Into<Vec<f32>>,
        closed: bool,
    ) -> Result<PiecewiseCubicCurve<f32>, PiecewiseMonotoneError> {
        let values = values.into();
        let slopes = vec![None; values.len()];
        PiecewiseCubicCurve::new_piecewise_monotone_with_slopes(values, slopes, grid, closed)
            .map_err(|e| match e {
                PiecewiseMonotoneWithSlopesError::FromPiecewiseMonotoneError(e) => e,
                _ => unreachable!(),
            })
    }
}

impl PiecewiseCubicCurve<f32> {
    pub fn new_piecewise_monotone_with_slopes(
        values: impl Into<Vec<f32>>,
        optional_slopes: impl AsRef<[Option<f32>]>,
        grid: impl Into<Vec<f32>>,
        closed: bool,
    ) -> Result<PiecewiseCubicCurve<f32>, PiecewiseMonotoneWithSlopesError> {
        use PiecewiseMonotoneError::*;
        use PiecewiseMonotoneWithSlopesError::*;
        let mut values = values.into();
        let optional_slopes = optional_slopes.as_ref();
        let mut grid = grid.into();
        if values.len() < 2 {
            return Err(LessThanTwoValues.into());
        }
        if values.len() != optional_slopes.len() {
            return Err(SlopesVsValues {
                slopes: optional_slopes.len(),
                values: values.len(),
            });
        }
        if values.len() + closed as usize != grid.len() {
            return Err(GridVsValues {
                grid: grid.len(),
                values: values.len(),
                closed,
            }
            .into());
        }
        check_grid(&grid).map_err(FromGridError)?;
        if closed {
            values.reserve_exact(2);
            grid.reserve_exact(1);
            // Closing the curve:
            values.push(values[0]);
            // One additional (temporary) value to calculate slopes:
            values.push(values[1]);
            grid.push(*grid.last().unwrap() + grid[1] - grid[0]);
        }
        let mut slopes = Vec::new();
        for i in 0..values.len() - 2 {
            if let (&[x_1, x0, x1, ..], &[t_1, t0, t1, ..]) = (&values[i..], &grid[i..]) {
                let left = (x0 - x_1) / (t0 - t_1);
                let right = (x1 - x0) / (t1 - t0);
                let slope = match optional_slopes[(i + 1) % optional_slopes.len()] {
                    Some(slope) => verify_slope(slope, left, right, i + 1)?,
                    None => fix_slope(catmull_rom_slope([x_1, x0, x1], [t_1, t0, t1]), left, right),
                };
                slopes.push(slope); // incoming
                slopes.push(slope); // outgoing
            } else {
                unreachable!();
            }
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
            slopes.push(calculate_slope(one, two, chord, 0)?);
            slopes.push(calculate_slope(two, one, chord, 1)?);
        } else {
            slopes.insert(
                0,
                calculate_slope(
                    *optional_slopes.first().unwrap(),
                    slopes.first().copied(),
                    (values[1] - values[0]) / (grid[1] - grid[0]),
                    0,
                )?,
            );
            slopes.push(calculate_slope(
                *optional_slopes.last().unwrap(),
                slopes.last().copied(),
                (values[values.len() - 1] - values[values.len() - 2])
                    / (grid[grid.len() - 1] - grid[grid.len() - 2]),
                slopes.len(),
            )?);
        }
        PiecewiseCubicCurve::new_hermite(&values, &slopes, &grid).map_err(|e| {
            use crate::cubichermitespline::Error as E;
            match e {
                E::LessThanTwoPositions => unreachable!(),
                E::TangentsVsSegments { .. } => unreachable!(),
                E::GridVsPositions { .. } => unreachable!(),
                E::FromGridError(e) => FromGridError(e).into(),
            }
        })
    }
}

pub(crate) fn catmull_rom_slope([x_1, x0, x1]: [f32; 3], [t_1, t0, t1]: [f32; 3]) -> f32 {
    let v0 = (x1 - x0) / (t1 - t0);
    let v_1 = (x0 - x_1) / (t0 - t_1);
    ((t1 - t0) * v_1 + (t0 - t_1) * v0) / (t1 - t_1)
}

fn calculate_slope(
    main: Option<f32>,
    other: Option<f32>,
    chord: f32,
    index: usize,
) -> Result<f32, PiecewiseMonotoneWithSlopesError> {
    if let Some(main) = main {
        Ok(verify_slope(main, chord, chord, index)?)
    } else if let Some(other) = other {
        Ok(end_slope(other, chord))
    } else {
        Ok(chord)
    }
}

fn verify_slope(
    slope: f32,
    left: f32,
    right: f32,
    index: usize,
) -> Result<f32, PiecewiseMonotoneWithSlopesError> {
    use PiecewiseMonotoneWithSlopesError::*;
    let fixed_slope = fix_slope(slope, left, right);
    #[allow(clippy::float_cmp)]
    if fixed_slope == slope {
        Ok(slope)
    } else if left * right < 0.0 {
        assert!(fixed_slope == 0.0);
        Err(SlopeTooSteep {
            index,
            slope,
            maximum: 0.0,
        })
    } else if left * slope < 0.0 {
        Err(SlopeWrongSign { index, slope })
    } else {
        Err(SlopeTooSteep {
            index,
            slope,
            maximum: fixed_slope,
        })
    }
}

/// Manipulate the slope to preserve shape.
/// See Dougherty et al. (1989), eq. (4.2).
pub(crate) fn fix_slope(slope: f32, left: f32, right: f32) -> f32 {
    if left * right <= 0.0 {
        0.0
    } else if right > 0.0 {
        slope.max(0.0).min(3.0 * left.abs().min(right.abs()))
    } else {
        slope.min(0.0).max(-3.0 * left.abs().min(right.abs()))
    }
}

/// NB: This is a very ad-hoc algorithm meant to minimize the change in slope
/// within the first/last curve segment.  Especially, this should avoid a
/// change from negative to positive acceleration (and vice versa).
/// There might be a better method available!?!
fn end_slope(inner_slope: f32, chord_slope: f32) -> f32 {
    if chord_slope < 0.0 {
        return -end_slope(-inner_slope, -chord_slope);
    }
    assert!(0.0 <= inner_slope);
    assert!(inner_slope <= 3.0 * chord_slope);
    if inner_slope <= chord_slope {
        3.0 * chord_slope - 2.0 * inner_slope
    } else {
        (3.0 * chord_slope - inner_slope) / 2.0
    }
}
