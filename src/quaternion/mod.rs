use std::borrow::Cow;

// Rename to avoid cbindgen error "'UnitQuaternion is not generic"
use nalgebra::UnitQuaternion as GenericUnitQuaternion;
use nalgebra::Vector3;

pub type UnitQuaternion = GenericUnitQuaternion<f32>;

pub type Vec3 = Vector3<f32>;

pub mod centripetalkochanekbartelsspline;
pub mod cubicdecasteljau;

pub use cubicdecasteljau::CubicDeCasteljau;

use crate::NormWrapper;

pub struct AngularVelocityNorm;

impl NormWrapper<AngularVelocityNorm> for Vec3 {
    fn norm(&self) -> f32 {
        self.norm()
    }
}

/// This is probably more complicated than it needs to be.
/// It's mainly for experimenting with the `Cow` type.
pub fn canonicalize<'a>(
    quaternions: impl Into<Cow<'a, [UnitQuaternion]>>,
) -> Cow<'a, [UnitQuaternion]> {
    let quaternions = quaternions.into();
    let mut p = UnitQuaternion::identity();

    let mut predicate = move |q: &UnitQuaternion| {
        let result = p.dot(q) < 0.0;
        p = *q;
        result
    };

    if let Some(idx) = quaternions.iter().position(predicate) {
        let mut quaternions = quaternions.into_owned();
        quaternions[idx].inverse_mut();
        for q in &mut quaternions[idx + 1..] {
            if predicate(q) {
                q.inverse_mut();
            }
        }
        quaternions.into()
    } else {
        quaternions
    }
}

/// angles in degrees!
pub fn angles2quat(azim: f32, elev: f32, roll: f32) -> UnitQuaternion {
    UnitQuaternion::from_axis_angle(&Vec3::z_axis(), azim.to_radians())
        * UnitQuaternion::from_axis_angle(&Vec3::x_axis(), elev.to_radians())
        * UnitQuaternion::from_axis_angle(&Vec3::y_axis(), roll.to_radians())
}
