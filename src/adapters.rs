use std::marker::PhantomData;

use crate::utilities::{bisect, GridError};
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
    #[error("first grid value must be specified")]
    FirstGridMissing,
    #[error("last grid value must be specified")]
    LastGridMissing,
    #[error("index {index}: duplicate value without corresponding grid value")]
    DuplicateValueWithoutGrid { index: usize },
    #[error(transparent)]
    FromGridError(#[from] GridError),
}

#[derive(thiserror::Error, Debug)]
pub enum NewGridWithSpeedsError {
    #[error(transparent)]
    FromNewGridError(#[from] NewGridError),
    #[error("number of grid values ({grid}) must be {} speeds ({speeds})", if *.closed {
        "one more than"
    } else {
        "the same as"
    })]
    GridVsSpeeds {
        grid: usize,
        speeds: usize,
        closed: bool,
    },
    #[error("index {index}: speed is only allowed if corresponding grid value is given")]
    SpeedWithoutGrid { index: usize },
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
        new_grid: impl AsRef<[Option<f32>]>,
        closed: bool,
    ) -> Result<NewGridAdapter<Value, Inner>, NewGridError> {
        let new_grid = new_grid.as_ref();
        let speeds = vec![None; inner.grid().len() - closed as usize];
        Self::adapt_with_speeds(inner, new_grid, speeds, closed).map_err(|e| match e {
            NewGridWithSpeedsError::FromNewGridError(e) => e,
            _ => unreachable!(),
        })
    }

    pub fn adapt_with_speeds(
        inner: Inner,
        new_grid: impl AsRef<[Option<f32>]>,
        speeds: impl AsRef<[Option<f32>]>,
        closed: bool,
    ) -> Result<NewGridAdapter<Value, Inner>, NewGridWithSpeedsError> {
        use NewGridError::*;
        use NewGridWithSpeedsError::*;
        let new_grid = new_grid.as_ref();
        let speeds = speeds.as_ref();
        if new_grid.len() != inner.grid().len() {
            return Err(NewGridVsOldGrid {
                new: new_grid.len(),
                old: inner.grid().len(),
            }
            .into());
        }
        if new_grid.len() != speeds.len() + closed as usize {
            return Err(GridVsSpeeds {
                grid: new_grid.len(),
                speeds: speeds.len(),
                closed,
            });
        }

        let mut t2u_times = Vec::new();
        let mut t2u_speeds = Vec::new();
        let mut missing_times = Vec::new();
        if let Some(time) = new_grid[0] {
            t2u_times.push(time);
        } else {
            return Err(FirstGridMissing.into());
        }
        t2u_speeds.push(speeds[0]);
        for i in 1..speeds.len() {
            let speed = speeds[i];
            if let Some(time) = new_grid[i] {
                t2u_times.push(time);
                t2u_speeds.push(speed);
            } else if speed.is_none() {
                missing_times.push(i);
            } else {
                return Err(SpeedWithoutGrid { index: i });
            }
        }
        if let Some(last_time) = *new_grid.last().unwrap() {
            if closed {
                t2u_times.push(last_time);
                t2u_speeds.push(speeds[0]);
            } else {
                // The last values have already been pushed in the for-loop above.
            }
        } else {
            return Err(LastGridMissing.into());
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
        let cyclic = closed && t2u_speeds[0].is_none();
        if cyclic {
            assert!(matches!(t2u_speeds[..], [None, .., None]));
        }
        let t2u = MonotoneCubicSpline::with_slopes(u_grid, t2u_speeds, t2u_times, cyclic).map_err(
            |e| {
                let fix_index = |mut idx| {
                    for &i in &missing_times {
                        if idx >= i {
                            idx += 1;
                        } else {
                            break;
                        }
                    }
                    idx
                };

                use crate::monotonecubicspline::MonotoneWithSlopesError as E;
                match e {
                    E::FromMonotoneError(e) => {
                        use crate::monotonecubicspline::MonotoneError as E;
                        match e {
                            // TODO: this might actually happen?
                            E::LessThanTwoValues => unreachable!(),
                            E::GridVsValues { .. } => unreachable!(),
                            E::Decreasing => unreachable!(),
                            E::FromGridError(mut e) => {
                                use crate::utilities::GridError::*;
                                match e {
                                    GridNan { ref mut index } => *index = fix_index(*index),
                                    GridNotAscending { ref mut index } => {
                                        *index = fix_index(*index)
                                    }
                                };
                                NewGridError::from(e).into()
                            }
                        }
                    }
                    E::SlopesVsValues { .. } => unreachable!(),
                    E::CyclicWithSlope { .. } => unreachable!(),
                    E::SlopeTooSteep {
                        index,
                        slope,
                        maximum,
                    } => TooFast {
                        index: fix_index(index),
                        speed: slope,
                        maximum,
                    },
                    E::NegativeSlope { index, slope } => NegativeSpeed {
                        index: fix_index(index),
                        speed: slope,
                    },
                }
            },
        )?;
        assert!(missing_times.len() == u_missing.len());
        for i in 0..missing_times.len() {
            if let Some(time) = t2u.get_time(u_missing[i]) {
                grid.insert(missing_times[i], time);
            } else {
                return Err(DuplicateValueWithoutGrid { index: i }.into());
            }
        }
        let t2u = t2u.into_inner();
        assert_eq!(inner.grid().len(), grid.len());

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
