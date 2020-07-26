use num_traits::one;

use crate::utilities::{check_grid, GridError};
use crate::{Scalar, Spline, SplineWithVelocity, Vector};

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

    pub fn segments(&self) -> &[[V; 4]] {
        &self.segments
    }

    // If t is out of bounds, it is trimmed to the smallest/largest possible value
    fn get_segment(&self, t: S) -> (S, S, S, &[V; 4]) {
        let (t, idx) = self.clamp_parameter_and_find_index(t);
        (t, self.grid[idx], self.grid[idx + 1], &self.segments[idx])
    }
}

impl<S: Scalar, V: Vector<S>> Spline<S, V> for PiecewiseCubicCurve<S, V> {
    fn evaluate(&self, t: S) -> V {
        let (t, t0, t1, a) = self.get_segment(t);
        let t = (t - t0) / (t1 - t0);
        ((a[3] * t + a[2]) * t + a[1]) * t + a[0]
    }

    fn grid(&self) -> &[S] {
        &self.grid
    }
}

impl<S, V> SplineWithVelocity<S, V, V> for PiecewiseCubicCurve<S, V>
where
    S: Scalar,
    V: Vector<S>,
{
    fn evaluate_velocity(&self, t: S) -> V {
        let (t, t0, t1, a) = self.get_segment(t);
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

    use crate::NormWrapper;

    struct NormF32;

    impl NormWrapper<NormF32> for f32 {
        type Norm = f32;

        fn norm(&self) -> Self::Norm {
            self.abs()
        }
    }

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
        assert_eq!(curve.integrated_speed::<NormF32>(0, 5.0, 6.0), 9.5);
        assert_eq!(curve.integrated_speed::<NormF32>(0, 5.0, 5.5), 2.5);
    }

    #[test]
    #[should_panic(expected = "assertion failed")]
    fn segment_length_early_begin() {
        let curve = make_simple_curve();
        curve.integrated_speed::<NormF32>(0, 4.9, 5.5);
    }

    #[test]
    #[should_panic(expected = "assertion failed")]
    fn segment_length_late_end() {
        let curve = make_simple_curve();
        curve.integrated_speed::<NormF32>(0, 5.1, 6.1);
    }

    #[test]
    fn grid() {
        let curve = make_simple_curve();
        assert_eq!(curve.grid(), &[5.0, 6.0]);
    }
}
