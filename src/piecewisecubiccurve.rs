use crate::utilities::{check_grid, GridError};
use crate::{Spline, SplineWithVelocity, Vector};

pub struct PiecewiseCubicCurve<V> {
    segments: Box<[[V; 4]]>,
    grid: Box<[f32]>,
}

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("there must be at least one segment")]
    ZeroSegments,
    #[error("length of grid ({grid}) must be one more than number of segments ({segments})")]
    GridVsSegments { grid: usize, segments: usize },
    #[error(transparent)]
    FromGridError(#[from] GridError),
}

impl<V: Vector> PiecewiseCubicCurve<V> {
    pub fn new(
        segments: impl Into<Box<[[V; 4]]>>,
        grid: impl Into<Box<[f32]>>,
    ) -> Result<PiecewiseCubicCurve<V>, Error> {
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
    fn get_segment(&self, t: f32) -> (f32, f32, f32, &[V; 4]) {
        let (t, idx) = self.clamp_parameter_and_find_index(t);
        (t, self.grid[idx], self.grid[idx + 1], &self.segments[idx])
    }
}

impl<V: Vector> Spline<V> for PiecewiseCubicCurve<V> {
    fn evaluate(&self, t: f32) -> V {
        let (t, t0, t1, a) = self.get_segment(t);
        let t = (t - t0) / (t1 - t0);
        ((a[3] * t + a[2]) * t + a[1]) * t + a[0]
    }

    fn grid(&self) -> &[f32] {
        &self.grid
    }
}

impl<V> SplineWithVelocity<V, V> for PiecewiseCubicCurve<V>
where
    V: Vector,
{
    fn evaluate_velocity(&self, t: f32) -> V {
        let (t, t0, t1, a) = self.get_segment(t);
        let t = (t - t0) / (t1 - t0);
        ((a[3] * 3.0 * t + a[2] * 2.0) * t + a[1]) / (t1 - t0)
    }
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    use crate::NormWrapper;

    struct NormF32;

    impl NormWrapper<NormF32> for f32 {
        fn norm(&self) -> f32 {
            self.abs()
        }
    }

    fn make_simple_curve() -> PiecewiseCubicCurve<f32> {
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
