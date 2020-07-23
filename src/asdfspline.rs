use crate::adapters::{ConstantSpeedAdapter, NewGridAdapter};
use crate::{PiecewiseCubicCurve, Scalar, Vector};

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

impl<S: Scalar> From<crate::adapters::NewGridWithSpeedsError<S>> for Error<S> {
    fn from(err: crate::adapters::NewGridWithSpeedsError<S>) -> Self {
        use crate::adapters::NewGridWithSpeedsError as Other;
        use Error::*;
        match err {
            Other::FirstTimeMissing => FirstTimeMissing,
            Other::LastTimeMissing => LastTimeMissing,
            Other::SpeedWithoutTime { index } => SpeedWithoutTime { index },
            Other::DuplicatePositionWithoutTime { index } => DuplicatePositionWithoutTime { index },
            Other::TimeNan { index } => TimeNan { index },
            // TODO: fix index?
            Other::GridNotAscending { index } => TimesNotAscending { index },
            // TODO: fix index?
            Other::TooFast {
                index,
                speed,
                maximum,
            } => TooFast {
                index,
                speed,
                maximum,
            },
            // TODO: fix index?
            Other::NegativeSpeed { index, speed } => NegativeSpeed { index, speed },
        }
    }
}

pub type AsdfPosSpline<S, V, F> =
    NewGridAdapter<S, V, ConstantSpeedAdapter<S, V, V, PiecewiseCubicCurve<S, V>, F>>;

impl<S, V, F> AsdfPosSpline<S, V, F>
where
    S: Scalar,
    V: Vector<S>,
    F: Fn(V) -> S,
{
    pub fn new(
        positions: &[V],
        times: &[Option<S>],
        speeds: &[Option<S>],
        tcb: &[[S; 3]],
        closed: bool,
        get_length: F,
    ) -> Result<AsdfPosSpline<S, V, F>, Error<S>> {
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
        let constant_speed = ConstantSpeedAdapter::adapt(path, get_length);
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

    type AsdfPosSpline1<F> = AsdfPosSpline<f32, f32, F>;

    fn get_length(x: f32) -> f32 {
        x.abs()
    }

    #[test]
    fn simple_linear() {
        let s = AsdfPosSpline1::new(
            &[1.0, 2.0],
            &[Some(0.0), Some(3.0)],
            &[None, None],
            &[],
            false,
            get_length,
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
            get_length,
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
            get_length,
        )
        .unwrap();
        assert_eq!(s.evaluate(4.0), 2.0);
    }
}
