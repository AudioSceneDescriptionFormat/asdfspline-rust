use num_traits::zero;
use superslice::Ext; // for slice::upper_bound_by()

use crate::utilities::bisect;
use crate::{MonotoneCubicSpline, PiecewiseCubicCurve, Scalar, Vector};

#[derive(thiserror::Error, Debug)]
pub enum Error<S: Scalar> {
    #[error("there must be at least two positions")]
    LessThanTwoPositions,
    #[error("number of times ({times}) must be {} number of positions ({positions})", if *.closed {
        "one more than"
    } else {
        "the same as"
    })]
    TimesVsPositions {
        times: usize,
        positions: usize,
        closed: bool,
    },
    #[error("number of speeds ({speeds}) and positions ({positions}) must be the same")]
    SpeedsVsPositions { speeds: usize, positions: usize },
    #[error("index {index}: speed is only allowed if time is given")]
    SpeedWithoutTime { index: usize },
    #[error("last time value must be specified")]
    LastTimeMissing,
    // TODO: two indices?
    #[error("index {index}: duplicate position without time")]
    DuplicatePositionWithoutTime { index: usize },
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
    #[error("index {index}: time values are not allowed to be NaN")]
    TimeNan { index: usize },
    #[error("index {index}: time values must be strictly ascending")]
    TimesNotAscending { index: usize },
    #[error("speed at index {index} too fast ({speed:?}; maximum: {maximum:?})")]
    TooFast { index: usize, speed: S, maximum: S },
    #[error("negative speed ({speed:?}) at index {index}")]
    NegativeSpeed { index: usize, speed: S },
}

impl<S: Scalar> From<crate::centripetalkochanekbartelsspline::Error> for Error<S> {
    fn from(err: crate::centripetalkochanekbartelsspline::Error) -> Self {
        use crate::centripetalkochanekbartelsspline::Error as Other;
        use Error::*;
        match err {
            Other::LessThanTwoPositions => LessThanTwoPositions,
            // TODO: fix index? consider missing times?
            // TODO: move to map_err() for this?
            // TODO: can there be non-adjacent indices?
            Other::RepeatedPosition { index } => RepeatedPosition { index },
            // TODO: same questions as above
            Other::TcbVsPositions {
                tcb,
                positions,
                closed,
            } => TcbVsPositions {
                tcb,
                positions,
                closed,
            },
        }
    }
}

pub struct AsdfSpline<S, V> {
    path: PiecewiseCubicCurve<S, V>,
    t2s: PiecewiseCubicCurve<S, S>, // Created via MonotoneCubicSpline
    grid: Box<[S]>,
    s_grid: Box<[S]>,
}

impl<S: Scalar, V: Vector<S>> AsdfSpline<S, V> {
    pub fn new<F: Fn(V) -> S>(
        positions: &[V],
        times: &[Option<S>],
        speeds: &[Option<S>],
        tcb: &[[S; 3]],
        closed: bool,
        get_length: F,
    ) -> Result<AsdfSpline<S, V>, Error<S>> {
        use Error::*;
        if positions.len() + closed as usize != times.len() {
            return Err(TimesVsPositions {
                times: times.len(),
                positions: positions.len(),
                closed,
            });
        }
        if speeds.len() != positions.len() {
            return Err(SpeedsVsPositions {
                speeds: speeds.len(),
                positions: positions.len(),
            });
        }
        let path = PiecewiseCubicCurve::new_centripetal_kochanek_bartels(
            positions,
            tcb,
            closed,
            &get_length,
        )?;

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
                return Err(SpeedWithoutTime { index: i });
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
            return Err(LastTimeMissing);
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
        use crate::monotonecubicspline::Error as Other;
        let t2s = MonotoneCubicSpline::with_slopes(lengths, t2s_speeds, t2s_times).map_err(
            |e| match e {
                Other::LessThanTwoValues => unreachable!(),
                Other::SlopesVsValues { .. } => unreachable!(),
                Other::GridVsValues { .. } => unreachable!(),
                Other::Decreasing => unreachable!(),
                // TODO: fix index? consider missing times? optional speeds?
                Other::GridNan { index } => TimeNan { index },
                // TODO: fix index?
                Other::GridNotAscending { index } => TimesNotAscending { index },
                // TODO: fix index?
                Other::SlopeTooSteep {
                    index,
                    slope,
                    maximum,
                } => TooFast {
                    index,
                    speed: slope,
                    maximum,
                },
                // TODO: fix index?
                Other::NegativeSlope { index, slope } => NegativeSpeed {
                    index,
                    speed: slope,
                },
            },
        )?;
        assert!(missing_times.len() == lengths_at_missing_times.len());
        for i in 0..missing_times.len() {
            if let Some(time) = t2s.get_time(lengths_at_missing_times[i]) {
                grid.insert(missing_times[i], time);
            } else {
                return Err(DuplicatePositionWithoutTime { index: i });
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
            &[1.0, 2.0],
            &[None, Some(3.0)],
            &[None, None],
            &[],
            false,
            get_length,
        )
        .unwrap();
        assert_eq!(s.evaluate(1.5, get_length), 1.5);
    }

    #[test]
    fn simple_closed() {
        let s = AsdfSpline1::new(
            &[1.0, 2.0],
            &[None, None, Some(3.0)],
            &[None, None],
            &[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            true,
            get_length,
        )
        .unwrap();
        assert_eq!(s.evaluate(1.5, get_length), 2.0);
    }

    #[test]
    fn closed_with_time() {
        let s = AsdfSpline1::new(
            &[1.0, 2.0],
            &[Some(3.0), Some(4.0), Some(5.0)],
            &[None, None],
            &[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            true,
            get_length,
        )
        .unwrap();
        assert_eq!(s.evaluate(4.0, get_length), 2.0);
    }
}
