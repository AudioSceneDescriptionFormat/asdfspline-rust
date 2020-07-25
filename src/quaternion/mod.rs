// Rename to avoid cbindgen error "'UnitQuaternion is not generic"
use nalgebra::UnitQuaternion as GenericUnitQuaternion;
use nalgebra::Vector3;

type UnitQuaternion = GenericUnitQuaternion<f32>;

type Vec3 = Vector3<f32>;

pub mod cubicdecasteljau;
