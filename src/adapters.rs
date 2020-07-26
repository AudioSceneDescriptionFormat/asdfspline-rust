use std::marker::PhantomData;

use num_traits::zero;

use crate::utilities::bisect;
use crate::{
    MonotoneCubicSpline, NormWrapper, PiecewiseCubicCurve, Scalar, Spline, SplineWithVelocity,
    Vector,
};

pub struct ConstantSpeedAdapter<S, Output, Velocity, Inner, U>
where
    S: Scalar,
    Velocity: Vector<S> + NormWrapper<U, Norm = S>,
    Inner: SplineWithVelocity<S, Output, Velocity>,
{
    inner: Inner,
    grid: Box<[S]>,
    _phantom_output: PhantomData<Output>,
    _phantom_velocity: PhantomData<Velocity>,
    _phantom_dummy: PhantomData<U>,
}

impl<S, Output, Velocity, Inner, U> ConstantSpeedAdapter<S, Output, Velocity, Inner, U>
where
    S: Scalar,
    Velocity: Vector<S> + NormWrapper<U, Norm = S>,
    Inner: SplineWithVelocity<S, Output, Velocity>,
{
    pub fn adapt(inner: Inner) -> ConstantSpeedAdapter<S, Output, Velocity, Inner, U> {
        let mut grid = Vec::with_capacity(inner.grid().len());
        grid.push(zero());
        let grid = inner
            .grid()
            .windows(2)
            .enumerate()
            .fold(grid, |mut l, (i, ts)| {
                if let [t0, t1] = *ts {
                    l.push(*l.last().unwrap() + inner.integrated_speed(i, t0, t1));
                    l
                } else {
                    unreachable!()
                }
            });
        ConstantSpeedAdapter {
            inner,
            grid: grid.into(),
            _phantom_output: PhantomData,
            _phantom_velocity: PhantomData,
            _phantom_dummy: PhantomData,
        }
    }

    /// If s is outside, return clipped t.
    /// This is only supposed to be used with `f32`.
    fn s2t(&self, s: S) -> S {
        // TODO: proper accuracy (a bit less than single-precision?)
        // TODO: a separate version for f64?
        let accuracy = S::from_f32(0.0001).unwrap();
        let (s, idx) = self.clamp_parameter_and_find_index(s);
        let mut s = s;
        s -= self.grid[idx];
        let t0 = self.inner.grid()[idx];
        let t1 = self.inner.grid()[idx + 1];
        let func = |t| self.inner.integrated_speed(idx, t0, t) - s;
        bisect(func, t0, t1, accuracy, 50)
    }
}

impl<S, Output, Velocity, Inner, U> Spline<S, Output>
    for ConstantSpeedAdapter<S, Output, Velocity, Inner, U>
where
    S: Scalar,
    Velocity: Vector<S> + NormWrapper<U, Norm = S>,
    Inner: SplineWithVelocity<S, Output, Velocity>,
{
    fn evaluate(&self, s: S) -> Output {
        self.inner.evaluate(self.s2t(s))
    }

    fn grid(&self) -> &[S] {
        &self.grid
    }
}

pub struct NewGridAdapter<S, Output, Inner>
where
    S: Scalar,
    Inner: Spline<S, Output>,
{
    inner: Inner,
    grid: Box<[S]>,
    t2u: PiecewiseCubicCurve<S, S>, // Created via MonotoneCubicSpline
    _phantom_output: PhantomData<Output>,
}

#[derive(thiserror::Error, Debug)]
pub enum NewGridWithSpeedsError<S: Scalar> {
    #[error("first time value must be specified")]
    FirstTimeMissing,
    #[error("last time value must be specified")]
    LastTimeMissing,
    #[error("index {index}: speed is only allowed if time is given")]
    SpeedWithoutTime { index: usize },
    #[error("index {index}: duplicate position without time")]
    DuplicatePositionWithoutTime { index: usize },
    #[error("index {index}: time values are not allowed to be NaN")]
    TimeNan { index: usize },
    #[error("index {index}: grid values must be strictly ascending")]
    GridNotAscending { index: usize },
    #[error("speed at index {index} too fast ({speed:?}; maximum: {maximum:?})")]
    TooFast { index: usize, speed: S, maximum: S },
    #[error("negative speed ({speed:?}) at index {index}")]
    NegativeSpeed { index: usize, speed: S },
}

impl<S, Output, Inner> NewGridAdapter<S, Output, Inner>
where
    S: Scalar,
    Inner: Spline<S, Output>,
{
    pub fn adapt_with_speeds(
        inner: Inner,
        times: &[Option<S>],
        speeds: &[Option<S>],
        closed: bool,
    ) -> Result<NewGridAdapter<S, Output, Inner>, NewGridWithSpeedsError<S>> {
        use NewGridWithSpeedsError::*;

        // TODO: check length of times and speeds, compare with inner.grid, check "closed"

        let mut t2u_times = Vec::new();
        let mut t2u_speeds = Vec::new();
        let mut missing_times = Vec::new();
        if let Some(time) = times[0] {
            t2u_times.push(time);
        } else {
            return Err(FirstTimeMissing);
        }
        t2u_speeds.push(speeds[0]);
        for i in 1..speeds.len() {
            let speed = speeds[i];
            if let Some(time) = times[i] {
                t2u_times.push(time);
                t2u_speeds.push(speed);
            } else if speed.is_none() {
                missing_times.push(i);
            } else {
                return Err(SpeedWithoutTime { index: i });
            }
        }
        if let Some(last_time) = *times.last().unwrap() {
            if closed {
                assert!(speeds.len() + 1 == times.len());
                t2u_times.push(last_time);
                t2u_speeds.push(speeds[0]);
            } else {
                // The last values have already been pushed in the for-loop above.
            }
        } else {
            return Err(LastTimeMissing);
        }

        let mut u_grid = Vec::<S>::new();
        let mut u_missing = Vec::<S>::new();
        for (i, &u) in inner.grid().iter().enumerate() {
            if missing_times.iter().any(|&x| x == i) {
                u_missing.push(u);
            } else {
                u_grid.push(u);
            }
        }

        let mut grid = t2u_times.clone();
        use crate::monotonecubicspline::Error as Other;
        let t2u =
            MonotoneCubicSpline::with_slopes(u_grid, t2u_speeds, t2u_times).map_err(
                |e| match e {
                    Other::LessThanTwoValues => unreachable!(),
                    Other::SlopesVsValues { .. } => unreachable!(),
                    Other::GridVsValues { .. } => unreachable!(),
                    Other::Decreasing => unreachable!(),
                    // TODO: fix index? consider missing times? optional speeds?
                    Other::GridNan { index } => TimeNan { index },
                    // TODO: fix index?
                    Other::GridNotAscending { index } => GridNotAscending { index },
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
        assert!(missing_times.len() == u_missing.len());
        for i in 0..missing_times.len() {
            if let Some(time) = t2u.get_time(u_missing[i]) {
                grid.insert(missing_times[i], time);
            } else {
                return Err(DuplicatePositionWithoutTime { index: i });
            }
        }
        let t2u = t2u.into_inner();
        assert_eq!(inner.grid().len(), grid.len());

        //let s_grid = grid.iter().map(|&t| t2u.evaluate(t)).collect();

        Ok(NewGridAdapter {
            inner,
            t2u,
            grid: grid.into(),
            _phantom_output: PhantomData,
        })
    }
}

impl<S, Output, Inner> Spline<S, Output> for NewGridAdapter<S, Output, Inner>
where
    S: Scalar,
    Inner: Spline<S, Output>,
{
    fn evaluate(&self, t: S) -> Output {
        self.inner.evaluate(self.t2u.evaluate(t))
    }

    fn grid(&self) -> &[S] {
        &self.grid
    }
}

// TODO: implement SplineWithVelocity?
/*
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
*/
