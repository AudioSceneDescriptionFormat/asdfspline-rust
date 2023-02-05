use std::borrow::Cow;

use superslice::Ext; // for slice::equal_range_by()

use crate::utilities::{bisect, check_grid, GridError};
use crate::PiecewiseCubicCurve;
use crate::Spline;
use crate::SplineWithVelocity;

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("values must not be decreasing")]
    Decreasing,
    #[error("there must be at least two values")]
    LessThanTwoValues,
    #[error("number of slopes ({slopes}) must be same as values ({values})")]
    SlopesVsValues { slopes: usize, values: usize },
    #[error("length of grid ({grid}) must be the same as number of values ({values})")]
    GridVsValues { grid: usize, values: usize },
    #[error(transparent)]
    FromGridError(#[from] GridError),
    #[error("slope at index {index} too steep ({slope:?}; maximum: {maximum:?})")]
    SlopeTooSteep {
        index: usize,
        slope: f32,
        maximum: f32,
    },
    #[error("negative slope ({slope:?}) at index {index}")]
    NegativeSlope { index: usize, slope: f32 },
    #[error(
        "If 'cyclic', the first and last slope must be (None, None), not ({first:?}, {last:?})"
    )]
    CyclicWithSlope {
        first: Option<f32>,
        last: Option<f32>,
    },
}

impl From<crate::shapepreservingcubicspline::Error> for Error {
    fn from(e: crate::shapepreservingcubicspline::Error) -> Self {
        use crate::shapepreservingcubicspline::Error as Other;
        use Error::*;
        match e {
            Other::LessThanTwoValues => LessThanTwoValues,
            Other::SlopesVsValues { slopes, values } => SlopesVsValues { slopes, values },
            Other::GridVsValues {
                grid,
                values,
                closed: false,
            } => GridVsValues { grid, values },
            Other::GridVsValues { closed: true, .. } => unreachable!(),
            Other::FromGridError(e) => e.into(),
            Other::SlopeTooSteep {
                index,
                slope,
                maximum,
            } => SlopeTooSteep {
                index,
                slope,
                maximum,
            },
            Other::SlopeWrongSign { index, slope } => NegativeSlope { index, slope },
        }
    }
}

/// ... monotonically *increasing* ...
pub struct MonotoneCubicSpline {
    inner: PiecewiseCubicCurve<f32>,
    values: Box<[f32]>,
}

impl MonotoneCubicSpline {
    pub fn new(
        values: impl Into<Box<[f32]>>,
        grid: impl Into<Vec<f32>>,
        cyclic: bool,
    ) -> Result<MonotoneCubicSpline, Error> {
        let values = values.into();
        let slopes = vec![None; values.len()];
        MonotoneCubicSpline::with_slopes(values, slopes, grid, cyclic)
    }

    pub fn with_slopes<'a>(
        values: impl Into<Box<[f32]>>,
        optional_slopes: impl Into<Cow<'a, [Option<f32>]>>,
        grid: impl Into<Vec<f32>>,
        cyclic: bool,
    ) -> Result<MonotoneCubicSpline, Error> {
        use Error::*;
        let values = values.into();
        let mut optional_slopes = optional_slopes.into();
        let grid = grid.into();
        // TODO: use is_sorted() once it is stabilized
        //if !values.is_sorted() {
        if values.windows(2).any(|w| w[0] > w[1]) {
            return Err(Decreasing);
        }

        // TODO: code re-use with new_shape_preserving_with_slopes()?
        if values.len() < 2 {
            return Err(LessThanTwoValues);
        }
        if values.len() != optional_slopes.len() {
            return Err(SlopesVsValues {
                slopes: optional_slopes.len(),
                values: values.len(),
            });
        }
        if values.len() != grid.len() {
            return Err(GridVsValues {
                grid: grid.len(),
                values: values.len(),
            });
        }
        check_grid(&grid)?;

        let closed = false;
        if cyclic {
            match optional_slopes[..] {
                [None, .., None] => {}
                [first, .., last] => {
                    return Err(CyclicWithSlope { first, last });
                }
                [] | [_] => unreachable!(),
            }
            let temp_values = vec![
                values[values.len() - 2],
                values[values.len() - 1],
                values[values.len() - 1] + (values[1] - values[0]),
            ];
            let no_slopes = vec![None, None, None];
            let temp_grid = vec![
                grid[grid.len() - 2],
                grid[grid.len() - 1],
                grid[grid.len() - 1] + (grid[1] - grid[0]),
            ];
            let temp_spline = PiecewiseCubicCurve::new_shape_preserving_with_slopes(
                temp_values,
                no_slopes,
                temp_grid,
                closed,
            )?;
            let cyclic_slope = temp_spline.evaluate_velocity(temp_spline.grid()[1]);
            if let [ref mut first, .., ref mut last] = optional_slopes.to_mut()[..] {
                *first = Some(cyclic_slope);
                *last = Some(cyclic_slope);
            } else {
                unreachable!();
            }
        }
        Ok(MonotoneCubicSpline {
            inner: PiecewiseCubicCurve::new_shape_preserving_with_slopes(
                values.as_ref(), // NB: Values are copied
                optional_slopes,
                grid,
                closed,
            )?,
            values,
        })
    }

    #[must_use]
    pub fn inner_ref(&self) -> &PiecewiseCubicCurve<f32> {
        &self.inner
    }

    #[must_use]
    pub fn into_inner(self) -> PiecewiseCubicCurve<f32> {
        self.inner
    }

    /// Get the time instance for the given value.
    ///
    /// If the solution is not unique, `None` is returned.
    /// If "value" is outside the range, the first/last time is returned.
    // TODO: rename to something with "solve"?
    #[must_use]
    pub fn get_time(&self, value: f32) -> Option<f32> {
        // NB: If initially given values are monotone (which we checked above!),
        // repetitions (i.e. a plateau) can only occur at those exact values.

        // NB: This doesn't work if a value is NaN (but we've checked for that already
        let range = self
            .values
            .equal_range_by(|x| x.partial_cmp(&value).unwrap());

        if range.end == 0 {
            // Value too small
            Some(*self.inner.grid().first().unwrap())
        } else if range.start == self.values.len() {
            // Value too large
            Some(*self.inner.grid().last().unwrap())
        } else if range.len() == 1 {
            // Exactly one match
            Some(self.inner.grid()[range.start])
        } else if range.len() > 1 {
            // Multiple matches
            None
        } else {
            // Within range, but no exact match

            let idx = range.end - 1;
            let mut a = self.inner.segments()[idx];
            a[0] -= value;

            let time = bisect(
                |t| ((a[3] * t + a[2]) * t + a[1]) * t + a[0],
                0.0,
                1.0,
                // TODO: proper tolerance value
                0.0001,
                500,
            );
            assert!((0.0..=1.0).contains(&time));
            let t0 = self.inner.grid()[idx];
            let t1 = self.inner.grid()[idx + 1];
            Some(time * (t1 - t0) + t0)
        }
    }
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    #[test]
    fn get_time_linear() {
        let values = [1.0, 2.0].to_vec();
        let grid = [3.0, 4.0].to_vec();
        let cyclic = false;
        let spline = MonotoneCubicSpline::new(values, grid, cyclic).unwrap();
        assert_eq!(spline.get_time(0.5).unwrap(), 3.0);
        assert_eq!(spline.get_time(1.0).unwrap(), 3.0);
        assert_eq!(spline.get_time(1.5).unwrap(), 3.5);
        assert_eq!(spline.get_time(2.0).unwrap(), 4.0);
        assert_eq!(spline.get_time(2.5).unwrap(), 4.0);
    }

    #[test]
    fn get_time_plateau() {
        let values = [1.0, 2.0, 2.0, 3.0].to_vec();
        let grid = [4.0, 5.0, 6.0, 7.0].to_vec();
        let cyclic = false;
        let spline = MonotoneCubicSpline::new(values, grid, cyclic).unwrap();
        assert_eq!(spline.get_time(2.0), None);
    }
}
