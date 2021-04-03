use crate::{PiecewiseCubicCurve, Vector};

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
    pub fn new_centripetal_kochanek_bartels<F>(
        positions: &[V],
        tcb: &[[f32; 3]],
        closed: bool,
        norm: F,
    ) -> Result<PiecewiseCubicCurve<V>, Error>
    where
        F: Fn(&V) -> f32,
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
            let delta = norm(&(x1 - x0)).sqrt();
            if delta == 0.0 {
                return Err(RepeatedPosition { index: i + 1 });
            }
            grid.push(*grid.last().unwrap() + delta);
        }
        let mut tangents = Vec::<V>::new();
        assert_eq!(positions.len(), grid.len());
        assert_eq!(positions.len(), tcb.len() + 2);
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

            let denominator = (t1 - t0) * (t0 - t_1) * (t1 - t_1);
            let incoming = ((x0 - x_1) * c * (t1 - t0).powi(2)
                + (x1 - x0) * d * (t0 - t_1).powi(2))
                / denominator;
            let outgoing = ((x0 - x_1) * a * (t1 - t0).powi(2)
                + (x1 - x0) * b * (t0 - t_1).powi(2))
                / denominator;
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
            assert_eq!(grid.len(), 2);
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

            if let (&[x0, x1, ..], &[t0, t1, ..]) = (positions, &grid[..]) {
                tangents.insert(0, natural_end_tangent(x0, x1, t0, t1, tangents[0]));
            } else {
                unreachable!();
            }
            if let (&[.., x0, x1], &[.., t0, t1]) = (positions, &grid[..]) {
                tangents.push(natural_end_tangent(
                    x0,
                    x1,
                    t0,
                    t1,
                    *tangents.last().unwrap(),
                ));
            } else {
                unreachable!();
            }
        }
        use crate::cubichermitespline::Error as Other;
        PiecewiseCubicCurve::new_hermite(&positions, &tangents, &grid).map_err(|e| match e {
            Other::LessThanTwoPositions => unreachable!(),
            Other::TangentsVsSegments { .. } => unreachable!(),
            Other::GridVsPositions { .. } => unreachable!(),
            Other::FromGridError(..) => unreachable!(),
        })
    }
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    use crate::Spline; // for evaluate(), grid()

    #[test]
    fn test_1d() {
        let positions = [1.0f32, 2.0, 3.0].to_vec();
        let tcb = [[4.0, 5.0, 6.0]];
        let closed = false;
        let curve =
            PiecewiseCubicCurve::new_centripetal_kochanek_bartels(&positions, &tcb, closed, |x| {
                x.abs()
            })
            .unwrap();
        assert_eq!(curve.grid()[0], 0.0);
        assert_eq!(curve.evaluate(0.0), 1.0);
        assert_eq!(curve.evaluate(*curve.grid().last().unwrap()), 3.0);
    }
}
