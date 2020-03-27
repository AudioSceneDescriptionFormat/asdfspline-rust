use num_traits::{one, zero};
use superslice::Ext; // for slice::equal_range_by()

use crate::utilities::bisect;
use crate::{fail, Error, Scalar};
use crate::{NewShapePreservingWithSlopes, PiecewiseCubicCurve};

/// ... monotonically *increasing* ...
pub struct MonotoneCubicSpline<S> {
    inner: PiecewiseCubicCurve<S, S>,
    values: Box<[S]>,
}

impl<S: Scalar> MonotoneCubicSpline<S> {
    pub fn new(
        values: impl Into<Box<[S]>>,
        grid: impl Into<Vec<S>>,
    ) -> Result<MonotoneCubicSpline<S>, Error> {
        let values = values.into();
        let slopes = vec![None; values.len()];
        MonotoneCubicSpline::with_slopes(values, slopes, grid)
    }

    pub fn with_slopes(
        values: impl Into<Box<[S]>>,
        optional_slopes: impl AsRef<[Option<S>]>,
        grid: impl Into<Vec<S>>,
    ) -> Result<MonotoneCubicSpline<S>, Error> {
        let values = values.into();
        // TODO: us is_sorted() once it is stabilized
        //if !values.is_sorted() {
        if values.windows(2).any(|w| w[0] > w[1]) {
            fail!("Values must be increasing");
        }
        let closed = false;
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

    pub fn inner_ref(&self) -> &PiecewiseCubicCurve<S, S> {
        &self.inner
    }

    pub fn into_inner(self) -> PiecewiseCubicCurve<S, S> {
        self.inner
    }

    /// Get the time instance for the given value.
    ///
    /// If the solution is not unique, `None` is returned.
    /// If "value" is outside the range, the first/last time is returned.
    // TODO: rename to something with "solve"?
    pub fn get_time(&self, value: S) -> Option<S> {
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
                zero(),
                one(),
                // TODO: proper tolerance value
                S::from_f64(0.0001).unwrap(),
                500,
            );
            assert!(zero::<S>() <= time && time <= one());
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
        let spline = MonotoneCubicSpline::<f32>::new(values, grid).unwrap();
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
        let spline = MonotoneCubicSpline::<f32>::new(values, grid).unwrap();
        assert_eq!(spline.get_time(2.0), None);
    }
}
