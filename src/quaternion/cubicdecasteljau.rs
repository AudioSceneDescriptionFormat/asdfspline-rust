use crate::utilities::{check_grid, GridError};
use crate::{Spline, SplineWithVelocity};

use super::{UnitQuaternion, Vec3};

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("there must be at least two grid elements")]
    GridTooShort,
    #[error("there must be a control quaternion for each grid point, \
        plus two between each pair of grid points \
        ({grid} + 2 * ({grid} - 1) = {} != {control_polygon})",
        .grid + 2 * (.grid - 1))]
    GridVsControlPolygon { grid: usize, control_polygon: usize },
    #[error("index {index}: NaN values are not allowed in grid")]
    GridNan { index: usize },
    #[error("index {index}: grid values must be strictly ascending")]
    GridNotAscending { index: usize },
}

impl From<GridError> for Error {
    fn from(e: GridError) -> Error {
        use Error::*;
        match e {
            GridError::GridNan { index } => GridNan { index },
            GridError::GridNotAscending { index } => GridNotAscending { index },
        }
    }
}

pub struct CubicDeCasteljau {
    control_polygon: Box<[UnitQuaternion]>,
    grid: Box<[f32]>,
}

impl CubicDeCasteljau {
    pub fn new(
        control_polygon: impl Into<Box<[UnitQuaternion]>>,
        grid: impl Into<Box<[f32]>>,
    ) -> Result<CubicDeCasteljau, Error> {
        let control_polygon = control_polygon.into();
        let grid = grid.into();
        use Error::*;
        if grid.len() < 2 {
            return Err(GridTooShort);
        }
        if grid.len() + 2 * (grid.len() - 1) != control_polygon.len() {
            return Err(GridVsControlPolygon {
                grid: grid.len(),
                control_polygon: control_polygon.len(),
            });
        }
        check_grid(&grid)?;
        Ok(CubicDeCasteljau {
            control_polygon,
            grid,
        })
    }

    /// Applies two levels of Slerp until only two quaternions are left
    #[allow(clippy::many_single_char_names)]
    fn partial_de_casteljau(&self, t: f32) -> (UnitQuaternion, UnitQuaternion, f32, f32) {
        let (t, idx) = self.clamp_parameter_and_find_index(t);
        let t0 = self.grid[idx];
        let t1 = self.grid[idx + 1];
        let delta_t = t1 - t0;
        let t = (t - t0) / delta_t;

        let a = &self.control_polygon[idx * 3];
        let b = &self.control_polygon[idx * 3 + 1];
        let c = &self.control_polygon[idx * 3 + 2];
        let d = &self.control_polygon[idx * 3 + 3];

        // NB: slerp() panics if angle is 180 degrees!

        let ab = a.slerp(b, t);
        let bc = b.slerp(c, t);
        let cd = c.slerp(d, t);

        (ab.slerp(&bc, t), bc.slerp(&cd, t), t, delta_t)
    }
}

impl Spline<UnitQuaternion> for CubicDeCasteljau {
    fn evaluate(&self, t: f32) -> UnitQuaternion {
        let (one, two, t, _) = self.partial_de_casteljau(t);
        one.slerp(&two, t)
    }

    fn grid(&self) -> &[f32] {
        &self.grid
    }
}

impl SplineWithVelocity<UnitQuaternion, Vec3> for CubicDeCasteljau {
    fn evaluate_velocity(&self, t: f32) -> Vec3 {
        let (one, two, _, delta_t) = self.partial_de_casteljau(t);
        const DEGREE: f32 = 3.0; // cubic
        one.rotation_to(&two).scaled_axis() * DEGREE / delta_t
    }
}
