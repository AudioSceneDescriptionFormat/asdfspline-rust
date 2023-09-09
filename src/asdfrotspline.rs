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
        match err {
            Other::LessThanTwoQuaternions => Self::LessThanTwoQuaternions,
            Other::TcbVsQuaternions {
                tcb,
                quaternions,
                closed,
            } => Self::TcbVsQuaternions {
                tcb,
                quaternions,
                closed,
            },
            Other::RepeatedQuaternion { index } => Self::RepeatedQuaternion { index },
        }
    }
}

impl From<crate::utilities::GridError> for Error {
    fn from(err: crate::utilities::GridError) -> Self {
        use crate::utilities::GridError as Other;
        match err {
            Other::GridNan { index } => Self::TimeNan { index },
            Other::GridNotAscending { index } => Self::TimesNotAscending { index },
        }
    }
}

impl From<crate::adapters::NewGridError> for Error {
    fn from(err: crate::adapters::NewGridError) -> Self {
        use crate::adapters::NewGridError as Other;
        match err {
            Other::FirstTimeMissing => Self::FirstTimeMissing,
            Other::LastTimeMissing => Self::LastTimeMissing,
            Other::DuplicateValueWithoutTime { index } => {
                Self::DuplicateQuaternionWithoutTime { index }
            }
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
        times: impl AsRef<[Option<f32>]>,
        tcb: impl AsRef<[[f32; 3]]>,
        closed: bool,
    ) -> Result<AsdfRotSpline, Error> {
        let quaternions = quaternions.into();
        let times = times.as_ref();
        let tcb = tcb.as_ref();
        if quaternions.len() + closed as usize != times.len() {
            return Err(Error::TimesVsQuaternions {
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
