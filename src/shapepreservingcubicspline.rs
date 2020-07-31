use crate::utilities::{check_grid, GridError};
use crate::PiecewiseCubicCurve;

// TODO: Two error types? One doesn't need "slopes" errors.
#[derive(thiserror::Error, Debug)]
pub enum Error {
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
    #[error("index {index}: NaN values are not allowed in grid")]
    GridNan { index: usize },
    #[error("index {index}: grid values must be strictly ascending")]
    GridNotAscending { index: usize },
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

impl From<GridError> for Error {
    fn from(e: GridError) -> Error {
        use Error::*;
        match e {
            GridError::GridNan { index } => GridNan { index },
            GridError::GridNotAscending { index } => GridNotAscending { index },
        }
    }
}

impl PiecewiseCubicCurve<f32> {
    pub fn new_shape_preserving(
        values: impl Into<Vec<f32>>,
        grid: impl Into<Vec<f32>>,
        closed: bool,
    ) -> Result<PiecewiseCubicCurve<f32>, Error> {
        let values = values.into();
        let slopes = vec![None; values.len()];
        PiecewiseCubicCurve::new_shape_preserving_with_slopes(values, slopes, grid, closed)
    }
}

impl PiecewiseCubicCurve<f32> {
    pub fn new_shape_preserving_with_slopes(
        values: impl Into<Vec<f32>>,
        optional_slopes: impl AsRef<[Option<f32>]>,
        grid: impl Into<Vec<f32>>,
        closed: bool,
    ) -> Result<PiecewiseCubicCurve<f32>, Error> {
        use Error::*;
        let mut values = values.into();
        let optional_slopes = optional_slopes.as_ref();
        let mut grid = grid.into();
        if values.len() < 2 {
            return Err(LessThanTwoValues);
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
            });
        }
        check_grid(&grid)?;
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
            let x_1 = values[i];
            let x0 = values[i + 1];
            let x1 = values[i + 2];
            let t_1 = grid[i];
            let t0 = grid[i + 1];
            let t1 = grid[i + 2];
            let left = (x0 - x_1) / (t0 - t_1);
            let right = (x1 - x0) / (t1 - t0);
            let slope = match optional_slopes[(i + 1) % optional_slopes.len()] {
                Some(slope) => verify_slope(slope, left, right, i + 1)?,
                None => fix_slope(
                    ((x0 - x_1) / (t0 - t_1) + (x1 - x0) / (t1 - t0)) / 2.0,
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
            slopes.push(calculate_slope(&one, &two, chord, 0)?);
            slopes.push(calculate_slope(&two, &one, chord, 1)?);
        } else {
            slopes.insert(
                0,
                calculate_slope(
                    optional_slopes.first().unwrap(),
                    &slopes.first().copied(),
                    (values[1] - values[0]) / (grid[1] - grid[0]),
                    0,
                )?,
            );
            slopes.push(calculate_slope(
                optional_slopes.last().unwrap(),
                &slopes.last().copied(),
                (values[values.len() - 1] - values[values.len() - 2])
                    / (grid[grid.len() - 1] - grid[grid.len() - 2]),
                slopes.len(),
            )?);
        }
        use crate::cubichermitespline::Error as Other;
        PiecewiseCubicCurve::new_hermite(&values, &slopes, &grid).map_err(|e| match e {
            Other::LessThanTwoPositions => unreachable!(),
            Other::TangentsVsSegments { .. } => unreachable!(),
            Other::GridVsPositions { .. } => unreachable!(),
            Other::GridNotAscending { index } => GridNotAscending { index },
            Other::GridNan { index } => GridNan { index },
        })
    }
}

fn calculate_slope(
    main: &Option<f32>,
    other: &Option<f32>,
    chord: f32,
    index: usize,
) -> Result<f32, Error> {
    if let Some(main) = main {
        Ok(verify_slope(*main, chord, chord, index)?)
    } else if let Some(other) = other {
        Ok(end_slope(*other, chord))
    } else {
        Ok(chord)
    }
}

fn verify_slope(slope: f32, left: f32, right: f32, index: usize) -> Result<f32, Error> {
    let fixed_slope = fix_slope(slope, left, right);
    #[allow(clippy::float_cmp)]
    if fixed_slope == slope {
        Ok(slope)
    } else if left * right < 0.0 {
        assert!(fixed_slope == 0.0);
        Err(Error::SlopeTooSteep {
            index,
            slope,
            maximum: 0.0,
        })
    } else if left * slope < 0.0 {
        Err(Error::SlopeWrongSign { index, slope })
    } else {
        Err(Error::SlopeTooSteep {
            index,
            slope,
            maximum: fixed_slope,
        })
    }
}

/// Manipulate the slope to preserve shape.
/// See Dougherty et al. (1989), eq. (4.2).
fn fix_slope(slope: f32, left: f32, right: f32) -> f32 {
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
