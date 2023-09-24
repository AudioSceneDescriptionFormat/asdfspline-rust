use super::{canonicalize, negate, CubicDeCasteljau, UnitQuaternion};

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("there must be at least two quaternions")]
    LessThanTwoQuaternions,
    #[error("number of quaternions ({quaternions}) must be {} TCB values ({tcb})", if *.closed {
        "the same as"
    } else {
        "two more than"
    })]
    TcbVsQuaternions {
        tcb: usize,
        quaternions: usize,
        closed: bool,
    },
    #[error("repeated quaternion (at index {index}) is not allowed")]
    RepeatedQuaternion { index: usize },
}

fn calculate_control_quaternions(
    qs: &[UnitQuaternion],
    ts: &[f32],
    tcb: &[f32],
) -> (UnitQuaternion, UnitQuaternion) {
    #[allow(non_snake_case)]
    if let ([q_1, q0, q1], [t_1, t0, t1], [T, C, B]) = (qs, ts, tcb) {
        let a = (1.0 - T) * (1.0 + C) * (1.0 + B);
        let b = (1.0 - T) * (1.0 - C) * (1.0 - B);
        let c = (1.0 - T) * (1.0 - C) * (1.0 + B);
        let d = (1.0 - T) * (1.0 + C) * (1.0 - B);

        let q_in = q_1.rotation_to(q0);
        let q_out = q0.rotation_to(q1);
        // w means omega (i.e. the angular velocity vector)
        let w_in = q_in.scaled_axis() / (t0 - t_1);
        let w_out = q_out.scaled_axis() / (t1 - t0);

        let w0 = |weight_in, weight_out| {
            (weight_in * (t1 - t0) * w_in + weight_out * (t0 - t_1) * w_out) / (t1 - t_1)
        };

        let degree = 3.0;
        (
            UnitQuaternion::from_scaled_axis(-w0(c, d) * (t0 - t_1) / degree) * q0,
            UnitQuaternion::from_scaled_axis(w0(a, b) * (t1 - t0) / degree) * q0,
        )
    } else {
        unreachable!();
    }
}

/// Calculates second control quaternion given first and third.
fn natural_control_quaternion(first: &UnitQuaternion, third: &UnitQuaternion) -> UnitQuaternion {
    first.rotation_to(third).powf(0.5) * first
}

impl CubicDeCasteljau {
    pub fn new_centripetal_kochanek_bartels(
        quaternions: impl Into<Vec<UnitQuaternion>>,
        tcb: &[[f32; 3]],
        closed: bool,
    ) -> Result<CubicDeCasteljau, Error> {
        use Error::*;
        let mut quaternions = quaternions.into();
        if quaternions.len() < 2 {
            return Err(LessThanTwoQuaternions);
        }
        if tcb.len() + 2 * !closed as usize != quaternions.len() {
            return Err(TcbVsQuaternions {
                tcb: tcb.len(),
                quaternions: quaternions.len(),
                closed,
            });
        }
        if closed {
            quaternions.push(quaternions[0]);
        }

        canonicalize(&mut quaternions);

        // Create grid with centripetal parameterization

        let mut grid = Vec::with_capacity(quaternions.len() + 2 * closed as usize);
        grid.push(0.0);
        for i in 0..quaternions.len() - 1 {
            if let [q0, q1] = &quaternions[i..i + 2] {
                let delta = q0.rotation_to(q1).angle().sqrt();
                if delta == 0.0 {
                    return Err(RepeatedQuaternion { index: i + 1 });
                }
                grid.push(*grid.last().unwrap() + delta);
            } else {
                unreachable!();
            }
        }

        if closed {
            if let (&[first, mut second, ..], &[.., mut penultimate, last]) =
                (&quaternions[..], &quaternions[..])
            {
                if penultimate.dot(&first) < 0.0 {
                    negate(&mut penultimate);
                }
                if last.dot(&second) < 0.0 {
                    negate(&mut second);
                }
                quaternions.insert(0, penultimate);
                quaternions.push(second);
            } else {
                unreachable!();
            }
            if let (&[first, second, ..], &[.., penultimate, last]) = (&grid[..], &grid[..]) {
                grid.insert(0, first - (last - penultimate));
                grid.push(last + (second - first));
            } else {
                unreachable!();
            }
        }

        let mut control_polygon = Vec::new();
        assert!(quaternions.len() == grid.len());
        for i in 0..quaternions.len() - 2 {
            let (incoming, outgoing) = calculate_control_quaternions(
                &quaternions[i..i + 3],
                &grid[i..i + 3],
                &tcb[i % tcb.len()],
            );
            control_polygon.push(incoming);
            control_polygon.push(quaternions[i + 1]);
            control_polygon.push(outgoing);
        }
        if closed {
            let _ = control_polygon.remove(0);
            let _ = control_polygon.pop();
            let _ = grid.remove(0);
            let _ = grid.pop();
        } else if control_polygon.is_empty() {
            // spherical linear interpolation between two quaternions
            assert_eq!(grid.len(), 2);
            assert!(tcb.is_empty());
            if let [q0, q1] = quaternions[..] {
                let offset = (q1 * q0.inverse()).powf(1.0 / 3.0); // "cubic" spline, degree 3
                control_polygon.push(q0);
                control_polygon.push(offset * q0);
                control_polygon.push(offset.inverse() * q1);
                control_polygon.push(q1);
            } else {
                unreachable!();
            }
        } else {
            if let ([first, ..], [third, ..]) = (&quaternions[..], &control_polygon[..]) {
                let second = natural_control_quaternion(first, third);
                control_polygon.insert(0, second);
                control_polygon.insert(0, *first);
            } else {
                unreachable!();
            }
            // Now counting from the end ...
            if let ([.., third], [.., first]) = (&control_polygon[..], &quaternions[..]) {
                let second = natural_control_quaternion(first, third);
                control_polygon.push(second);
                control_polygon.push(*first);
            } else {
                unreachable!();
            }
        }
        CubicDeCasteljau::new(control_polygon, grid).map_err(|e| {
            use super::cubicdecasteljau::Error as E;
            match e {
                E::GridTooShort => unreachable!(),
                E::GridVsControlPolygon { .. } => unreachable!(),
                E::FromGridError(_) => unreachable!(),
            }
        })
    }
}
