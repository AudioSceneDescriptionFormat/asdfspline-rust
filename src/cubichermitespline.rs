use num_traits::one;

use crate::PiecewiseCubicCurve;
use crate::{Scalar, Vector};

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("There must be at least two positions")]
    LessThanTwoPositions,
    #[error(
        "Exactly 2 tangents per segment are required \
            (got {segments} segments and {tangents} tangents)"
    )]
    TangentsVsSegments { tangents: usize, segments: usize },
    #[error("length of grid ({grid}) must be the same as number of positions ({positions})")]
    GridVsPositions { grid: usize, positions: usize },
    #[error("index {index}: NaN values are not allowed in grid")]
    GridNan { index: usize },
    #[error("index {index}: grid values must be strictly ascending")]
    GridNotAscending { index: usize },
}

impl<S: Scalar, V: Vector<S>> PiecewiseCubicCurve<S, V> {
    pub fn new_hermite(
        positions: &[V],
        tangents: &[V],
        grid: &[S],
    ) -> Result<PiecewiseCubicCurve<S, V>, Error> {
        use Error::*;
        if positions.len() < 2 {
            return Err(LessThanTwoPositions);
        }
        let segments_len = positions.len() - 1;
        if tangents.len() != 2 * segments_len {
            return Err(TangentsVsSegments {
                tangents: tangents.len(),
                segments: segments_len,
            });
        }
        if positions.len() != grid.len() {
            return Err(GridVsPositions {
                grid: grid.len(),
                positions: positions.len(),
            });
        }
        let mut segments = Vec::with_capacity(segments_len);
        for i in 0..segments_len {
            let x0 = positions[i];
            let x1 = positions[i + 1];
            let v0 = tangents[2 * i];
            let v1 = tangents[2 * i + 1];
            let t0 = grid[i];
            let t1 = grid[i + 1];
            let delta = t1 - t0;
            let one: S = one();
            let two = one + one;
            let three = two + one;

            // [a0]   [ 1,  0,          0,      0] [x0]
            // [a1] = [ 0,  0,      delta,      0] [x1]
            // [a2]   [-3,  3, -2 * delta, -delta] [v0]
            // [a3]   [ 2, -2,      delta,  delta] [v1]

            segments.push([
                x0,
                v0 * delta,
                x0 * -three + x1 * three - v0 * two * delta - v1 * delta,
                x0 * two - x1 * two + v0 * delta + v1 * delta,
            ]);
        }
        use crate::piecewisecubiccurve::Error as Other;
        PiecewiseCubicCurve::new(segments, grid).map_err(|err| match err {
            Other::ZeroSegments => unreachable!(),
            Other::GridVsSegments { .. } => unreachable!(),
            Other::GridNan { index } => GridNan { index },
            Other::GridNotAscending { index } => GridNotAscending { index },
        })
    }
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    #[test]
    fn test_1d() {
        let positions = [1.0, 2.0];
        let tangents = [3.0, 4.0];
        let grid = [5.0, 6.0];
        let curve = PiecewiseCubicCurve::new_hermite(&positions, &tangents, &grid).unwrap();
        assert_eq!(curve.grid(), &[5.0, 6.0]);
    }

    #[test]
    fn test_errors() {
        let result = PiecewiseCubicCurve::new_hermite(&[1.0], &[3.0, 4.0], &[5.0, 6.0]);
        assert!(result.is_err());
        let result = PiecewiseCubicCurve::new_hermite(&[1.0, 2.0], &[3.0], &[5.0, 6.0]);
        assert!(result.is_err());
        let result = PiecewiseCubicCurve::new_hermite(&[1.0, 2.0], &[3.0, 3.5, 4.0], &[5.0, 6.0]);
        assert!(result.is_err());
        let result = PiecewiseCubicCurve::new_hermite(&[1.0, 2.0], &[3.0, 4.0], &[5.0]);
        assert!(result.is_err());
    }
}
