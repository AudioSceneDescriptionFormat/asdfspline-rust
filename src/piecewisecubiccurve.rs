use num_traits::one;
use superslice::Ext; // for slice::upper_bound_by()

use crate::utilities::{check_grid, gauss_legendre13, GridError};
use crate::{Scalar, Vector};

pub struct PiecewiseCubicCurve<S, V> {
    segments: Box<[[V; 4]]>,
    grid: Box<[S]>,
}

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("there must be at least one segment")]
    ZeroSegments,
    #[error("length of grid ({grid}) must be one more than number of segments ({segments})")]
    GridVsSegments { grid: usize, segments: usize },
    #[error("index {index}: NaN values are not allowed in grid")]
    GridNan { index: usize },
    #[error("index {index}: grid values must be strictly ascending")]
    GridNotAscending { index: usize },
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

impl<S: Scalar, V: Vector<S>> PiecewiseCubicCurve<S, V> {
    pub fn new(
        segments: impl Into<Box<[[V; 4]]>>,
        grid: impl Into<Box<[S]>>,
    ) -> Result<PiecewiseCubicCurve<S, V>, Error> {
        let segments = segments.into();
        let grid = grid.into();
        use Error::*;
        if segments.len() < 1 {
            return Err(ZeroSegments);
        }
        if segments.len() + 1 != grid.len() {
            return Err(GridVsSegments {
                grid: grid.len(),
                segments: segments.len(),
            });
        }
        check_grid(&grid)?;
        Ok(PiecewiseCubicCurve { segments, grid })
    }

    pub fn evaluate(&self, t: S) -> V {
        let (t, t0, t1, a) = self.get_segment(t);
        let t = (t - t0) / (t1 - t0);
        ((a[3] * t + a[2]) * t + a[1]) * t + a[0]
    }

    pub fn evaluate_velocity(&self, t: S) -> V {
        let (t, t0, t1, coeffs) = self.get_segment(t);
        PiecewiseCubicCurve::segment_velocity(t, t0, t1, coeffs)
    }

    pub fn segments(&self) -> &[[V; 4]] {
        &self.segments
    }

    pub fn grid(&self) -> &[S] {
        &self.grid
    }

    pub fn segment_length<F>(&self, index: usize, get_length: F) -> S
    where
        F: Fn(V) -> S,
    {
        let t0 = self.grid[index];
        let t1 = self.grid[index + 1];
        self.segment_partial_length(index, t0, t1, get_length)
    }

    pub fn segment_partial_length<F>(&self, index: usize, a: S, b: S, get_length: F) -> S
    where
        F: Fn(V) -> S,
    {
        assert!(a <= b);
        let coeffs = self.segments[index];
        let t0 = self.grid[index];
        let t1 = self.grid[index + 1];
        assert!(t0 <= a);
        assert!(b <= t1);

        let speed = |t| get_length(PiecewiseCubicCurve::segment_velocity(t, t0, t1, &coeffs));
        gauss_legendre13(speed, a, b)
    }

    // If t is out of bounds, it is trimmed to the smallest/largest possible value
    fn get_segment(&self, mut t: S) -> (S, S, S, &[V; 4]) {
        let first = *self.grid.first().unwrap();
        let last = *self.grid.last().unwrap();
        let idx = if t < first {
            t = first;
            0
        } else if t < last {
            // NB: This doesn't work if a value is NaN (but we've checked for that in new())
            self.grid.upper_bound_by(|x| x.partial_cmp(&t).unwrap()) - 1
        } else if t == last {
            self.grid.len() - 2
        } else {
            t = last;
            self.segments.len() - 1
        };
        (t, self.grid[idx], self.grid[idx + 1], &self.segments[idx])
    }

    fn segment_velocity(t: S, t0: S, t1: S, a: &[V; 4]) -> V {
        let t = (t - t0) / (t1 - t0);
        let one: S = one();
        let two = one + one;
        let three = two + one;
        ((a[3] * three * t + a[2] * two) * t + a[1]) / (t1 - t0)
    }
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    fn make_simple_curve() -> PiecewiseCubicCurve<f32, f32> {
        PiecewiseCubicCurve {
            segments: Box::new([[1.0, 2.5, 3.0, 4.0]]),
            grid: Box::new([5.0, 6.0]),
        }
    }

    #[test]
    fn evaluate() {
        let curve = make_simple_curve();
        assert_eq!(curve.evaluate(4.5), 1.0); // t < first
        assert_eq!(curve.evaluate(5.0), 1.0); // t == first
        assert_eq!(curve.evaluate(5.5), 3.5); // t < last
        assert_eq!(curve.evaluate(6.0), 10.5); // t == last
        assert_eq!(curve.evaluate(6.5), 10.5); // last < t
    }

    #[test]
    fn evaluate_velocity() {
        let curve = make_simple_curve();
        assert_eq!(curve.evaluate_velocity(5.0), 2.5);
        assert_eq!(curve.evaluate_velocity(5.5), 8.5);
        assert_eq!(curve.evaluate_velocity(6.0), 20.5);
    }

    #[test]
    fn segment_length() {
        let curve = make_simple_curve();
        assert_eq!(curve.segment_length(0, |x| x.abs()), 9.5);
        assert_eq!(curve.segment_partial_length(0, 5.0, 5.5, |x| x.abs()), 2.5);
    }

    #[test]
    #[should_panic(expected = "assertion failed")]
    fn segment_length_early_begin() {
        let curve = make_simple_curve();
        curve.segment_partial_length(0, 4.9, 5.5, |x| x.abs());
    }

    #[test]
    #[should_panic(expected = "assertion failed")]
    fn segment_length_late_end() {
        let curve = make_simple_curve();
        curve.segment_partial_length(0, 5.1, 6.1, |x| x.abs());
    }

    #[test]
    fn grid() {
        let curve = make_simple_curve();
        assert_eq!(curve.grid(), &[5.0, 6.0]);
    }
}
