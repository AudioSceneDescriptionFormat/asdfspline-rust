use num_traits::{one, pow, zero};

use crate::make_cubic_hermite_spline;
use crate::PiecewiseCubicCurve;
use crate::{fail, Error};
use crate::{Scalar, Vector};

// TODO: require alga::linear::NormedSpace to get norm() method for vectors?
pub fn make_centripetal_kochanek_bartels_spline<S, V, F>(
    vertices: impl Into<Vec<V>>,
    tcb: impl AsRef<[[S; 3]]>,
    closed: bool,
    get_length: F,
) -> Result<PiecewiseCubicCurve<S, V>, Error>
where
    S: Scalar,
    V: Vector<S>,
    F: Fn(V) -> S,
{
    let mut vertices = vertices.into();
    if vertices.len() < 2 {
        fail!("At least two vertices are required");
    }
    if closed {
        vertices.push(vertices[0]);
        vertices.push(vertices[1]);
    }

    let tcb = tcb.as_ref();
    if tcb.len() + 2 != vertices.len() {
        fail!("There must be two more vertices than TCB values (except for closed curves)");
    }

    // Create grid with centripetal parametrization

    let mut grid = Vec::<S>::with_capacity(vertices.len());
    grid.push(zero());
    for i in 0..vertices.len() - 1 {
        let x0 = vertices[i];
        let x1 = vertices[i + 1];
        let delta = get_length(x1 - x0).sqrt();
        if delta == zero() {
            fail!("Repeated vertices are not possible");
        }
        grid.push(*grid.last().unwrap() + delta);
    }
    let mut tangents = Vec::<V>::new();
    assert!(vertices.len() == grid.len());
    assert!(vertices.len() == tcb.len() + 2);
    for i in 0..vertices.len() - 2 {
        let x_1 = vertices[i];
        let x0 = vertices[i + 1];
        let x1 = vertices[i + 2];
        let t_1 = grid[i];
        let t0 = grid[i + 1];
        let t1 = grid[i + 2];
        #[allow(non_snake_case)]
        let [T, C, B] = tcb[(i + closed as usize) % tcb.len()];
        let one = one::<S>();
        let a = (one - T) * (one + C) * (one + B);
        let b = (one - T) * (one - C) * (one - B);
        let c = (one - T) * (one - C) * (one + B);
        let d = (one - T) * (one + C) * (one - B);

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

        // Remove temporary vertex and grid elements:
        vertices.pop();
        grid.pop();
    } else if vertices.len() == 2 {
        // Straight line
        assert!(grid.len() == 2);
        assert!(tangents.is_empty());
        let tangent = (vertices[1] - vertices[0]) / (grid[1] - grid[0]);
        tangents.push(tangent);
        tangents.push(tangent);
    } else {
        // End conditions for non-closed curves
        assert!(tangents.len() >= 2);

        let one: S = one();
        let two = one + one;
        let three = two + one;

        // "natural" end conditions
        let natural_end_tangent = |x0, x1, t0, t1, inner_tangent| {
            let delta = t1 - t0;
            (x1 * three - x0 * three - inner_tangent * delta) / (two * delta)
        };

        tangents.insert(
            0,
            natural_end_tangent(vertices[0], vertices[1], grid[0], grid[1], tangents[0]),
        );
        let last = vertices.len() - 1;
        tangents.push(natural_end_tangent(
            vertices[last - 1],
            vertices[last],
            grid[last - 1],
            grid[last],
            *tangents.last().unwrap(),
        ));
    }
    make_cubic_hermite_spline(&vertices, &tangents, &grid)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_1d() {
        let vertices = [1.0f32, 2.0, 3.0].to_vec();
        let tcb = [[4.0, 5.0, 6.0]];
        let closed = false;
        let curve =
            make_centripetal_kochanek_bartels_spline(vertices, &tcb, closed, |x| x.abs()).unwrap();
        assert_eq!(curve.grid()[0], 0.0);
        assert_eq!(curve.evaluate(0.0), 1.0);
        assert_eq!(curve.evaluate(*curve.grid().last().unwrap()), 3.0);
    }
}
