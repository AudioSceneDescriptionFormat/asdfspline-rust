use num_traits::one;

use crate::PiecewiseCubicCurve;
use crate::{fail, Error};
use crate::{Scalar, Vector};

pub fn make_cubic_hermite_spline<S: Scalar, V: Vector<S>>(
    positions: &[V],
    tangents: &[V],
    grid: &[S],
) -> Result<PiecewiseCubicCurve<S, V>, Error> {
    if positions.len() < 2 {
        fail!("At least 2 positions are needed");
    }
    let segments_len = positions.len() - 1;
    if tangents.len() != 2 * segments_len {
        fail!("Exactly 2 tangents per segment are needed");
    }
    if positions.len() != grid.len() {
        fail!("As many grid times as positions are needed");
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
    PiecewiseCubicCurve::new(segments, grid)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_1d() {
        let positions = [1.0, 2.0];
        let tangents = [3.0, 4.0];
        let grid = [5.0, 6.0];
        let curve = make_cubic_hermite_spline(&positions, &tangents, &grid).unwrap();
        assert_eq!(curve.grid(), &[5.0, 6.0]);
    }

    #[test]
    fn test_errors() {
        let result = make_cubic_hermite_spline(&[1.0], &[3.0, 4.0], &[5.0, 6.0]);
        assert!(result.is_err());
        let result = make_cubic_hermite_spline(&[1.0, 2.0], &[3.0], &[5.0, 6.0]);
        assert!(result.is_err());
        let result = make_cubic_hermite_spline(&[1.0, 2.0], &[3.0, 3.5, 4.0], &[5.0, 6.0]);
        assert!(result.is_err());
        let result = make_cubic_hermite_spline(&[1.0, 2.0], &[3.0, 4.0], &[5.0]);
        assert!(result.is_err());
    }
}
