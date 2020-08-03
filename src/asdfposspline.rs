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

impl From<crate::centripetalkochanekbartelsspline::Error> for Error {
    fn from(err: crate::centripetalkochanekbartelsspline::Error) -> Self {
        use crate::centripetalkochanekbartelsspline::Error as Other;
        use Error::*;
        match err {
            Other::LessThanTwoPositions => LessThanTwoPositions,
            Other::RepeatedPosition { index } => RepeatedPosition { index },
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

impl From<crate::utilities::GridError> for Error {
    fn from(err: crate::utilities::GridError) -> Self {
        use crate::utilities::GridError as Other;
        use Error::*;
        match err {
            Other::GridNan { index } => TimeNan { index },
            Other::GridNotAscending { index } => TimesNotAscending { index },
        }
    }
}

impl From<crate::adapters::NewGridError> for Error {
    fn from(err: crate::adapters::NewGridError) -> Self {
        use crate::adapters::NewGridError as Other;
        use Error::*;
        match err {
            Other::FirstTimeMissing => FirstTimeMissing,
            Other::LastTimeMissing => LastTimeMissing,
            Other::DuplicateValueWithoutTime { index } => DuplicatePositionWithoutTime { index },
            Other::NewGridVsOldGrid { .. } => unreachable!(),
            Other::FromGridError(e) => e.into(),
        }
    }
}

impl From<crate::adapters::NewGridWithSpeedsError> for Error {
    fn from(err: crate::adapters::NewGridWithSpeedsError) -> Self {
        use crate::adapters::NewGridWithSpeedsError as Other;
        use Error::*;
        match err {
            Other::FromNewGridError(e) => e.into(),
            Other::SpeedWithoutTime { index } => SpeedWithoutTime { index },
            Other::TooFast {
                index,
                speed,
                maximum,
            } => TooFast {
                index,
                speed,
                maximum,
            },
            Other::NegativeSpeed { index, speed } => NegativeSpeed { index, speed },
            Other::GridVsSpeeds { .. } => unreachable!(),
        }
    }
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
        )?;
        let constant_speed = ConstantSpeedAdapter::adapt(path);
        Ok(NewGridAdapter::adapt_with_speeds(
            constant_speed,
            times,
            speeds,
            closed,
        )?)
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
        let s = AsdfPosSpline1::new(
            &[1.0, 2.0],
            &[Some(0.0), Some(3.0)],
            &[None, None],
            &[],
            false,
        )
        .unwrap();
        assert_eq!(s.evaluate(1.5), 1.5);
    }

    #[test]
    fn simple_closed() {
        let s = AsdfPosSpline1::new(
            &[1.0, 2.0],
            &[Some(0.0), None, Some(3.0)],
            &[None, None],
            &[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            true,
        )
        .unwrap();
        assert_eq!(s.evaluate(1.5), 2.0);
    }

    #[test]
    fn closed_with_time() {
        let s = AsdfPosSpline1::new(
            &[1.0, 2.0],
            &[Some(3.0), Some(4.0), Some(5.0)],
            &[None, None],
            &[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            true,
        )
        .unwrap();
        assert_eq!(s.evaluate(4.0), 2.0);
    }
}
