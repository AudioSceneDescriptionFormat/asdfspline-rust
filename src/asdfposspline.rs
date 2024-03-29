use crate::adapters::{ConstantSpeedAdapter, NewGridAdapter};
use crate::{NormWrapper, PiecewiseCubicCurve, Vector};

#[derive(thiserror::Error, Debug)]
pub enum Error {
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
    #[error("first time value must be specified")]
    FirstTimeMissing,
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
    TooFast {
        index: usize,
        speed: f32,
        maximum: f32,
    },
    #[error("negative speed ({speed:?}) at index {index}")]
    NegativeSpeed { index: usize, speed: f32 },
}

pub type AsdfPosSpline<V, U> =
    NewGridAdapter<V, ConstantSpeedAdapter<V, V, PiecewiseCubicCurve<V>, U>>;

impl<V, U> AsdfPosSpline<V, U>
where
    V: Vector + NormWrapper<U>,
{
    pub fn new(
        positions: impl AsRef<[V]>,
        times: impl AsRef<[Option<f32>]>,
        speeds: impl AsRef<[Option<f32>]>,
        tcb: impl AsRef<[[f32; 3]]>,
        closed: bool,
    ) -> Result<AsdfPosSpline<V, U>, Error> {
        use Error::*;
        let positions = positions.as_ref();
        let times = times.as_ref();
        let speeds = speeds.as_ref();
        let tcb = tcb.as_ref();
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
            NormWrapper::norm,
        )
        .map_err(|e| {
            use crate::centripetalkochanekbartelsspline::Error as E;
            match e {
                E::LessThanTwoPositions => LessThanTwoPositions,
                E::RepeatedPosition { index } => RepeatedPosition { index },
                E::TcbVsPositions {
                    tcb,
                    positions,
                    closed,
                } => TcbVsPositions {
                    tcb,
                    positions,
                    closed,
                },
            }
        })?;
        let constant_speed = ConstantSpeedAdapter::adapt(path);
        NewGridAdapter::adapt_with_speeds(constant_speed, times, speeds, closed).map_err(|e| {
            use crate::adapters::NewGridWithSpeedsError as E;
            match e {
                E::FromNewGridError(e) => {
                    use crate::adapters::NewGridError as E;
                    match e {
                        E::FirstGridMissing => FirstTimeMissing,
                        E::LastGridMissing => LastTimeMissing,
                        E::DuplicateValueWithoutGrid { index } => {
                            DuplicatePositionWithoutTime { index }
                        }
                        E::NewGridVsOldGrid { .. } => unreachable!(),
                        E::FromGridError(e) => {
                            use crate::utilities::GridError as E;
                            match e {
                                E::GridNan { index } => TimeNan { index },
                                E::GridNotAscending { index } => TimesNotAscending { index },
                            }
                        }
                    }
                }
                E::SpeedWithoutGrid { index } => SpeedWithoutTime { index },
                E::TooFast {
                    index,
                    speed,
                    maximum,
                } => TooFast {
                    index,
                    speed,
                    maximum,
                },
                E::NegativeSpeed { index, speed } => NegativeSpeed { index, speed },
                E::GridVsSpeeds { .. } => unreachable!(),
            }
        })
    }
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    use crate::Spline; // for evaluate()

    struct NormF32;

    impl NormWrapper<NormF32> for f32 {
        fn norm(&self) -> f32 {
            self.abs()
        }
    }

    type AsdfPosSpline1 = AsdfPosSpline<f32, NormF32>;

    #[test]
    fn simple_linear() {
        let s = AsdfPosSpline1::new([1.0, 2.0], [Some(0.0), Some(3.0)], [None, None], [], false)
            .unwrap();
        assert_eq!(s.evaluate(1.5), 1.5);
    }

    #[test]
    fn simple_closed() {
        let s = AsdfPosSpline1::new(
            [1.0, 2.0],
            [Some(0.0), None, Some(3.0)],
            [None, None],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            true,
        )
        .unwrap();
        assert_eq!(s.evaluate(1.5), 2.0);
    }

    #[test]
    fn closed_with_time() {
        let s = AsdfPosSpline1::new(
            [1.0, 2.0],
            [Some(3.0), Some(4.0), Some(5.0)],
            [None, None],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            true,
        )
        .unwrap();
        assert_eq!(s.evaluate(4.0), 2.0);
    }
}
