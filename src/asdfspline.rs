use num_traits::zero;
use superslice::Ext; // for slice::upper_bound_by()

use crate::utilities::bisect;
use crate::{
    fail, make_centripetal_kochanek_bartels_spline, Error, MonotoneCubicSpline,
    PiecewiseCubicCurve, Scalar, Vector,
};

pub struct AsdfSpline<S, V> {
    path: PiecewiseCubicCurve<S, V>,
    t2s: PiecewiseCubicCurve<S, S>, // Created via MonotoneCubicSpline
    grid: Box<[S]>,
    s_grid: Box<[S]>,
}

impl<S: Scalar, V: Vector<S>> AsdfSpline<S, V> {
    pub fn new<F: Fn(V) -> S>(
        positions: impl Into<Vec<V>>,
        times: impl AsRef<[Option<S>]>,
        speeds: impl AsRef<[Option<S>]>,
        tcb: impl AsRef<[[S; 3]]>,
        closed: bool,
        get_length: F,
    ) -> Result<AsdfSpline<S, V>, Error> {
        let positions = positions.into();
        if positions.len() < 2 {
            fail!("At least two positions are required");
        }
        let times = times.as_ref();
        if positions.len() + closed as usize != times.len() {
            fail!("Number of time values must be same as positions (one more for closed curves)");
        }
        let speeds = speeds.as_ref();
        if speeds.len() != positions.len() {
            fail!("Same number of speed values as positions are required");
        }
        let path = make_centripetal_kochanek_bartels_spline(positions, tcb, closed, &get_length)?;

        let mut t2s_times = Vec::new();
        let mut t2s_speeds = Vec::new();
        let mut missing_times = Vec::new();
        if let Some(time) = times[0] {
            t2s_times.push(time);
        } else {
            t2s_times.push(zero());
        }
        t2s_speeds.push(speeds[0]);
        for i in 1..speeds.len() {
            let speed = speeds[i];
            if let Some(time) = times[i] {
                t2s_times.push(time);
                t2s_speeds.push(speed);
            } else if speed.is_none() {
                missing_times.push(i);
            } else {
                fail!("Speed is only allowed if time is given");
            }
        }
        if let Some(last_time) = *times.last().unwrap() {
            if closed {
                assert!(speeds.len() + 1 == times.len());
                t2s_times.push(last_time);
                t2s_speeds.push(t2s_speeds[0]);
            } else {
                // The last values have already been pushed in the for-loop above.
            }
        } else {
            fail!("Last time value must be specified");
        }
        let mut lengths = Vec::<S>::new();
        let mut lengths_at_missing_times = Vec::<S>::new();
        lengths.push(zero());
        for i in 0..path.grid().len() - 1 {
            let length = path.segment_length(i, &get_length);
            if missing_times.iter().any(|&x| x == i) {
                lengths_at_missing_times.push(*lengths.last().unwrap());
                *lengths.last_mut().unwrap() += length;
            } else {
                lengths.push(*lengths.last().unwrap() + length);
            }
        }
        let mut grid = t2s_times.clone();
        let t2s = MonotoneCubicSpline::with_slopes(lengths, t2s_speeds, t2s_times)?;
        assert!(missing_times.len() == lengths_at_missing_times.len());
        for i in 0..missing_times.len() {
            if let Some(time) = t2s.get_time(lengths_at_missing_times[i]) {
                grid.insert(missing_times[i], time);
            } else {
                fail!("duplicate vertex without time");
            }
        }
        let t2s = t2s.into_inner();
        assert_eq!(path.grid().len(), grid.len());
        let s_grid = grid.iter().map(|&t| t2s.evaluate(t)).collect();
        Ok(AsdfSpline {
            path,
            t2s,
            grid: grid.into(),
            s_grid,
        })
    }

    pub fn evaluate<F>(&self, t: S, get_length: F) -> V
    where
        F: Fn(V) -> S,
    {
        self.path
            .evaluate(self.s2u(self.t2s.evaluate(t), get_length))
    }

    pub fn evaluate_velocity<F>(&self, t: S, get_length: F) -> V
    where
        F: Fn(V) -> S,
    {
        let speed = self.t2s.evaluate_velocity(t);
        let mut tangent = self
            .path
            .evaluate_velocity(self.s2u(self.t2s.evaluate(t), &get_length));
        let tangent_length = get_length(tangent);
        if tangent_length != zero() {
            tangent /= tangent_length;
        }
        tangent * speed
    }

    pub fn grid(&self) -> &[S] {
        &self.grid
    }

    /// If s is outside, return clipped u.
    /// This is only supposed to be used with `f32`.
    fn s2u<F>(&self, s: S, get_length: F) -> S
    where
        F: Fn(V) -> S,
    {
        // TODO: proper accuracy (a bit less than single-precision?)
        // TODO: a separate version for f64?
        let accuracy = S::from_f32(0.0001).unwrap();

        let index = if s <= *self.s_grid.first().unwrap() {
            return *self.path.grid().first().unwrap();
        } else if s < *self.s_grid.last().unwrap() {
            // NB: This doesn't work if a value is NaN (but this shouldn't happen anyway)
            self.s_grid.upper_bound_by(|x| x.partial_cmp(&s).unwrap()) - 1
        } else {
            return *self.path.grid().last().unwrap();
        };
        let mut s = s;
        s -= self.s_grid[index];
        let u0 = self.path.grid()[index];
        let u1 = self.path.grid()[index + 1];
        let func = |u| self.path.segment_partial_length(index, u0, u, &get_length) - s;
        bisect(func, u0, u1, accuracy, 50)
    }
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    type AsdfSpline1 = AsdfSpline<f32, f32>;

    fn get_length(x: f32) -> f32 {
        x.abs()
    }

    #[test]
    fn simple_linear() {
        let s = AsdfSpline1::new(
            [1.0, 2.0].to_vec(),
            [None, Some(3.0)],
            [None, None],
            [],
            false,
            get_length,
        )
        .unwrap();
        assert_eq!(s.evaluate(1.5, get_length), 1.5);
    }

    #[test]
    fn simple_closed() {
        let s = AsdfSpline1::new(
            [1.0, 2.0].to_vec(),
            [None, None, Some(3.0)],
            [None, None],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            true,
            get_length,
        )
        .unwrap();
        assert_eq!(s.evaluate(1.5, get_length), 2.0);
    }

    #[test]
    fn closed_with_time() {
        let s = AsdfSpline1::new(
            [1.0, 2.0].to_vec(),
            [Some(3.0), Some(4.0), Some(5.0)],
            [None, None],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            true,
            get_length,
        )
        .unwrap();
        assert_eq!(s.evaluate(4.0, get_length), 2.0);
    }
}
