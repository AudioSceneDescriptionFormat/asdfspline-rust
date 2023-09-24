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
        use Error::*;
        let quaternions = quaternions.into();
        let times = times.as_ref();
        let tcb = tcb.as_ref();
        if quaternions.len() + closed as usize != times.len() {
            return Err(TimesVsQuaternions {
                times: times.len(),
                quaternions: quaternions.len(),
                closed,
            });
        }
        let path = CubicDeCasteljau::new_centripetal_kochanek_bartels(quaternions, tcb, closed)
            .map_err(|e| {
                use crate::quaternion::centripetalkochanekbartelsspline::Error as E;
                match e {
                    E::LessThanTwoQuaternions => LessThanTwoQuaternions,
                    E::TcbVsQuaternions {
                        tcb,
                        quaternions,
                        closed,
                    } => TcbVsQuaternions {
                        tcb,
                        quaternions,
                        closed,
                    },
                    E::RepeatedQuaternion { index } => RepeatedQuaternion { index },
                }
            })?;
        let constant_speed = ConstantSpeedAdapter::adapt(path);
        NewGridAdapter::adapt(constant_speed, times, closed).map_err(|e| {
            use crate::adapters::NewGridError as E;
            match e {
                E::FirstGridMissing => FirstTimeMissing,
                E::LastGridMissing => LastTimeMissing,
                E::DuplicateValueWithoutGrid { index } => DuplicateQuaternionWithoutTime { index },
                E::FromGridError(e) => {
                    use crate::utilities::GridError as E;
                    match e {
                        E::GridNan { index } => TimeNan { index },
                        E::GridNotAscending { index } => TimesNotAscending { index },
                    }
                }
                E::NewGridVsOldGrid { .. } => unreachable!(),
            }
        })
    }
}
