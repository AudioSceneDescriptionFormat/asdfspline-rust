use std::marker::PhantomData;

use crate::utilities::bisect;
use crate::{
    MonotoneCubicSpline, NormWrapper, PiecewiseCubicCurve, Spline, SplineWithVelocity, Vector,
};

pub struct ConstantSpeedAdapter<Value, Velocity, Inner, U>
where
    Velocity: Vector + NormWrapper<U>,
    Inner: SplineWithVelocity<Value, Velocity>,
{
    inner: Inner,
    grid: Box<[f32]>,
    _phantom_output: PhantomData<Value>,
    _phantom_velocity: PhantomData<Velocity>,
    _phantom_dummy: PhantomData<U>,
}

impl<Value, Velocity, Inner, U> ConstantSpeedAdapter<Value, Velocity, Inner, U>
where
    Velocity: Vector + NormWrapper<U>,
    Inner: SplineWithVelocity<Value, Velocity>,
{
    pub fn adapt(inner: Inner) -> ConstantSpeedAdapter<Value, Velocity, Inner, U> {
        let mut grid = Vec::with_capacity(inner.grid().len());
        grid.push(0.0);
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
    fn s2t(&self, s: f32) -> f32 {
        // TODO: proper accuracy (a bit less than single-precision?)
        let accuracy = 0.0001;
        let (s, idx) = self.clamp_parameter_and_find_index(s);
        let mut s = s;
        s -= self.grid[idx];
        let t0 = self.inner.grid()[idx];
        let t1 = self.inner.grid()[idx + 1];
        let func = |t| self.inner.integrated_speed(idx, t0, t) - s;
        bisect(func, t0, t1, accuracy, 50)
    }
}

impl<Value, Velocity, Inner, U> Spline<Value> for ConstantSpeedAdapter<Value, Velocity, Inner, U>
where
    Velocity: Vector + NormWrapper<U>,
    Inner: SplineWithVelocity<Value, Velocity>,
{
    fn evaluate(&self, s: f32) -> Value {
        self.inner.evaluate(self.s2t(s))
    }

    fn grid(&self) -> &[f32] {
        &self.grid
    }
}

pub struct NewGridAdapter<Value, Inner>
where
    Inner: Spline<Value>,
{
    inner: Inner,
    grid: Box<[f32]>,
    t2u: PiecewiseCubicCurve<f32>, // Created via MonotoneCubicSpline
    _phantom_output: PhantomData<Value>,
}

#[derive(thiserror::Error, Debug)]
pub enum NewGridError {
    #[error("length of new grid ({new}) must be same as old grid ({old})")]
    NewGridVsOldGrid { new: usize, old: usize },
    #[error("first time value must be specified")]
    FirstTimeMissing,
    #[error("last time value must be specified")]
    LastTimeMissing,
    #[error("index {index}: duplicate value without time")]
    DuplicateValueWithoutTime { index: usize },
    #[error("index {index}: time values are not allowed to be NaN")]
    TimeNan { index: usize },
    #[error("index {index}: grid values must be strictly ascending")]
    GridNotAscending { index: usize },
}

#[derive(thiserror::Error, Debug)]
pub enum NewGridWithSpeedsError {
    #[error(transparent)]
    FromNewGridError(#[from] NewGridError),
    #[error("number of times ({times}) must be {} speeds ({speeds})", if *.closed {
        "one more than"
    } else {
        "the same as"
    })]
    TimesVsSpeeds {
        times: usize,
        speeds: usize,
        closed: bool,
    },
    #[error("index {index}: speed is only allowed if time is given")]
    SpeedWithoutTime { index: usize },
    #[error("speed at index {index} too fast ({speed:?}; maximum: {maximum:?})")]
    TooFast {
        index: usize,
        speed: f32,
        maximum: f32,
    },
    #[error("negative speed ({speed:?}) at index {index}")]
    NegativeSpeed { index: usize, speed: f32 },
}

impl<Value, Inner> NewGridAdapter<Value, Inner>
where
    Inner: Spline<Value>,
{
    pub fn adapt(
        inner: Inner,
        times: &[Option<f32>],
        closed: bool,
    ) -> Result<NewGridAdapter<Value, Inner>, NewGridError> {
        let speeds = vec![None; inner.grid().len() - closed as usize];
        Self::adapt_with_speeds(inner, times, speeds, closed).map_err(|e| match e {
            NewGridWithSpeedsError::FromNewGridError(e) => e,
            _ => unreachable!(),
        })
    }

    pub fn adapt_with_speeds(
        inner: Inner,
        times: impl AsRef<[Option<f32>]>,
        speeds: impl AsRef<[Option<f32>]>,
        closed: bool,
    ) -> Result<NewGridAdapter<Value, Inner>, NewGridWithSpeedsError> {
        use NewGridError::*;
        use NewGridWithSpeedsError::*;
        let times = times.as_ref();
        let speeds = speeds.as_ref();
        if times.len() != inner.grid().len() {
            return Err(NewGridVsOldGrid {
                new: times.len(),
                old: inner.grid().len(),
            }
            .into());
        }
        if times.len() != speeds.len() + closed as usize {
            return Err(TimesVsSpeeds {
                times: times.len(),
                speeds: speeds.len(),
                closed,
            });
        }

        let mut t2u_times = Vec::new();
        let mut t2u_speeds = Vec::new();
        let mut missing_times = Vec::new();
        if let Some(time) = times[0] {
            t2u_times.push(time);
        } else {
            return Err(FirstTimeMissing.into());
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
                t2u_times.push(last_time);
                t2u_speeds.push(speeds[0]);
            } else {
                // The last values have already been pushed in the for-loop above.
            }
        } else {
            return Err(LastTimeMissing.into());
        }

        let mut u_grid = Vec::new();
        let mut u_missing = Vec::new();
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
                    Other::GridNan { index } => TimeNan { index }.into(),
                    // TODO: fix index?
                    Other::GridNotAscending { index } => GridNotAscending { index }.into(),
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
                return Err(DuplicateValueWithoutTime { index: i }.into());
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

impl<Value, Inner> Spline<Value> for NewGridAdapter<Value, Inner>
where
    Inner: Spline<Value>,
{
    fn evaluate(&self, t: f32) -> Value {
        self.inner.evaluate(self.t2u.evaluate(t))
    }

    fn grid(&self) -> &[f32] {
        &self.grid
    }
}
