use num_traits::pow;

use crate::{NormWrapper, PiecewiseCubicCurve, Vector};

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("there must be at least two positions")]
    LessThanTwoPositions,
    #[error("number of positions ({positions}) must be {} TCB values ({tcb})", if *.closed {
        "the same as"
    } else {
        "two more than"
    })]
    TcbVsPositions {
        tcb: usize,
        positions: usize,
        closed: bool,
    },
    #[error("repeated position (at index {index}) is not allowed")]
    RepeatedPosition { index: usize },
}

impl<V: Vector> PiecewiseCubicCurve<V> {
    pub fn new_centripetal_kochanek_bartels<U>(
        positions: &[V],
        tcb: &[[f32; 3]],
        closed: bool,
    ) -> Result<PiecewiseCubicCurve<V>, Error>
    where
        V: NormWrapper<U>,
    {
        use Error::*;
        let positions_len = positions.len();
        if positions_len < 2 {
            return Err(LessThanTwoPositions);
        }
        let mut positions = positions;
        // Only used for "closed" splines:
        let mut positions_vec;
        if closed {
            positions_vec = Vec::with_capacity(positions.len() + 2);
            positions_vec.extend(positions);
            positions_vec.push(positions[0]);
            positions_vec.push(positions[1]);
            positions = &positions_vec;
        } else {
            // To avoid error: "use of possibly uninitialized `positions_vec`"
            positions_vec = Vec::new();
        }

        if tcb.len() + 2 != positions.len() {
            return Err(TcbVsPositions {
                tcb: tcb.len(),
                positions: positions_len,
                closed,
            });
        }

        // Create grid with centripetal parametrization

        let mut grid = Vec::with_capacity(positions.len());
        grid.push(0.0);
        for i in 0..positions.len() - 1 {
            let x0 = positions[i];
            let x1 = positions[i + 1];
            let delta = (x1 - x0).norm().sqrt();
            if delta == 0.0 {
                return Err(RepeatedPosition { index: i + 1 });
            }
            grid.push(*grid.last().unwrap() + delta);
        }
        let mut tangents = Vec::<V>::new();
        assert!(positions.len() == grid.len());
        assert!(positions.len() == tcb.len() + 2);
        for i in 0..positions.len() - 2 {
            let x_1 = positions[i];
            let x0 = positions[i + 1];
            let x1 = positions[i + 2];
            let t_1 = grid[i];
            let t0 = grid[i + 1];
            let t1 = grid[i + 2];
            #[allow(non_snake_case)]
            let [T, C, B] = tcb[(i + closed as usize) % tcb.len()];
            let a = (1.0 - T) * (1.0 + C) * (1.0 + B);
            let b = (1.0 - T) * (1.0 - C) * (1.0 - B);
            let c = (1.0 - T) * (1.0 - C) * (1.0 + B);
            let d = (1.0 - T) * (1.0 + C) * (1.0 - B);

            let incoming = ((x0 - x_1) * c * pow(t1 - t0, 2) + (x1 - x0) * d * pow(t0 - t_1, 2))
                / ((t1 - t0) * (t0 - t_1) * (t1 - t_1));
            let outgoing = ((x0 - x_1) * a * pow(t1 - t0, 2) + (x1 - x0) * b * pow(t0 - t_1, 2))
                / ((t1 - t0) * (t0 - t_1) * (t1 - t_1));
            tangents.push(incoming);
            tangents.push(outgoing);
        }

        if closed {
            // Move last (outgoing) tangent to the beginning:
            tangents.rotate_right(1);

            // Remove temporary position and grid elements:
            positions_vec.pop();
            grid.pop();

            // Update reference
            positions = &positions_vec;
        } else if positions.len() == 2 {
            // Straight line
            assert!(grid.len() == 2);
            assert!(tangents.is_empty());
            let tangent = (positions[1] - positions[0]) / (grid[1] - grid[0]);
            tangents.push(tangent);
            tangents.push(tangent);
        } else {
            // End conditions for non-closed curves
            assert!(tangents.len() >= 2);

            // "natural" end conditions
            let natural_end_tangent = |x0, x1, t0, t1, inner_tangent| {
                let delta = t1 - t0;
                (x1 * 3.0 - x0 * 3.0 - inner_tangent * delta) / (2.0 * delta)
            };

            tangents.insert(
                0,
                natural_end_tangent(positions[0], positions[1], grid[0], grid[1], tangents[0]),
            );
            let last = positions.len() - 1;
            tangents.push(natural_end_tangent(
                positions[last - 1],
                positions[last],
                grid[last - 1],
                grid[last],
                *tangents.last().unwrap(),
            ));
        }
        use crate::cubichermitespline::Error as Other;
        PiecewiseCubicCurve::new_hermite(&positions, &tangents, &grid).map_err(|e| match e {
            Other::LessThanTwoPositions => unreachable!(),
            Other::TangentsVsSegments { .. } => unreachable!(),
            Other::GridVsPositions { .. } => unreachable!(),
            Other::GridNotAscending { .. } => unreachable!(),
            Other::GridNan { .. } => unreachable!(),
        })
    }
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    use crate::Spline; // for evaluate(), grid()

    struct NormF32;

    impl NormWrapper<NormF32> for f32 {
        fn norm(&self) -> f32 {
            self.abs()
        }
    }

    #[test]
    fn test_1d() {
        let positions = [1.0f32, 2.0, 3.0].to_vec();
        let tcb = [[4.0, 5.0, 6.0]];
        let closed = false;
        let curve = PiecewiseCubicCurve::new_centripetal_kochanek_bartels::<NormF32>(
            &positions, &tcb, closed,
        )
        .unwrap();
        assert_eq!(curve.grid()[0], 0.0);
        assert_eq!(curve.evaluate(0.0), 1.0);
        assert_eq!(curve.evaluate(*curve.grid().last().unwrap()), 3.0);
    }
}
