use crate::adapters::{ConstantSpeedAdapter, NewGridAdapter};
use crate::quaternion::{AngularVelocityNorm, CubicDeCasteljau, UnitQuaternion, Vec3};

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("there must be at least two quaternions")]
    LessThanTwoQuaternions,
    #[error("number of times ({times}) must be {} number of quaternions ({quaternions})", if *.closed {
        "one more than"
    } else {
        "the same as"
    })]
    TimesVsQuaternions {
        times: usize,
        quaternions: usize,
        closed: bool,
    },
    #[error("first time value must be specified")]
    FirstTimeMissing,
    #[error("last time value must be specified")]
    LastTimeMissing,
    #[error("index {index}: duplicate quaternion without time")]
    DuplicateQuaternionWithoutTime { index: usize },
    #[error("index {index}: time values are not allowed to be NaN")]
    TimeNan { index: usize },
    #[error("index {index}: time values must be strictly ascending")]
    TimesNotAscending { index: usize },
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

impl From<crate::quaternion::centripetalkochanekbartelsspline::Error> for Error {
    fn from(err: crate::quaternion::centripetalkochanekbartelsspline::Error) -> Self {
        use crate::quaternion::centripetalkochanekbartelsspline::Error as Other;
        use Error::*;
        match err {
            Other::LessThanTwoQuaternions => LessThanTwoQuaternions,
            Other::TcbVsQuaternions {
                tcb,
                quaternions,
                closed,
            } => TcbVsQuaternions {
                tcb,
                quaternions,
                closed,
            },
            Other::RepeatedQuaternion { index } => RepeatedQuaternion { index },
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
            Other::DuplicateValueWithoutTime { index } => DuplicateQuaternionWithoutTime { index },
            Other::FromGridError(e) => e.into(),
            Other::NewGridVsOldGrid { .. } => unreachable!(),
        }
    }
}

pub type AsdfRotSpline = NewGridAdapter<
    UnitQuaternion,
    ConstantSpeedAdapter<UnitQuaternion, Vec3, CubicDeCasteljau, AngularVelocityNorm>,
>;

impl AsdfRotSpline {
    pub fn new(
        quaternions: impl Into<Vec<UnitQuaternion>>,
        times: &[Option<f32>],
        tcb: &[[f32; 3]],
        closed: bool,
    ) -> Result<AsdfRotSpline, Error> {
        use Error::*;
        let quaternions = quaternions.into();
        if quaternions.len() + closed as usize != times.len() {
            return Err(TimesVsQuaternions {
                times: times.len(),
                quaternions: quaternions.len(),
                closed,
            });
        }
        let path = CubicDeCasteljau::new_centripetal_kochanek_bartels(quaternions, tcb, closed)?;
        let constant_speed = ConstantSpeedAdapter::adapt(path);
        Ok(NewGridAdapter::adapt(constant_speed, times, closed)?)
    }
}
