// Re-export to make it easy for downstream crates to use the proper version:
pub use nalgebra;

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

pub fn canonicalize(quaternions: &mut [UnitQuaternion]) {
    let mut p = UnitQuaternion::identity();
    for q in quaternions {
        if p.dot(q) < 0.0 {
            q.inverse_mut();
        }
        p = *q;
    }
}

/// angles in degrees!
#[must_use]
pub fn angles2quat(azim: f32, elev: f32, roll: f32) -> UnitQuaternion {
    UnitQuaternion::from_axis_angle(&Vec3::z_axis(), azim.to_radians())
        * UnitQuaternion::from_axis_angle(&Vec3::x_axis(), elev.to_radians())
        * UnitQuaternion::from_axis_angle(&Vec3::y_axis(), roll.to_radians())
}
