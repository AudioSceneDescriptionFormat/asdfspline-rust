// Rename to avoid cbindgen error "'UnitQuaternion is not generic"
use nalgebra::UnitQuaternion as GenericUnitQuaternion;
use nalgebra::Vector3;

type UnitQuaternion = GenericUnitQuaternion<f32>;

type Vec3 = Vector3<f32>;

pub mod centripetalkochanekbartelsspline;
pub mod cubicdecasteljau;

pub use cubicdecasteljau::CubicDeCasteljau;

pub fn canonicalize(quaternions: &mut [UnitQuaternion]) {
    let mut p = UnitQuaternion::identity();
    for q in quaternions {
        if p.dot(&q) < 0.0 {
            p.inverse_mut();
        }
        p = *q;
    }
}
