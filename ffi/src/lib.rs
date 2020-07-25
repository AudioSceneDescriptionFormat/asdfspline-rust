use std::cell::RefCell;
use std::ffi::CString;
use std::fmt::Display;
use std::slice;

use libc::{c_char, size_t};
use nalgebra_glm as glm;

use asdfspline::{AsdfPosSpline, MonotoneCubicSpline, NormWrapper, PiecewiseCubicCurve, Spline};

thread_local! {
    static LAST_ERROR: RefCell<CString> = RefCell::new(CString::new("no error").unwrap());
}

fn set_error<D: Display>(error: D) {
    LAST_ERROR.with(|cell| {
        *cell.borrow_mut() = CString::new(error.to_string()).unwrap();
    });
}

/// The error message will be freed if another error occurs. It is the caller's
/// responsibility to make sure they're no longer using the string before
/// calling any other function which may fail.
#[no_mangle]
pub extern "C" fn asdf_last_error() -> *const c_char {
    LAST_ERROR.with(|cell| cell.borrow().as_ptr())
}

trait ResultExt<T, E: Display> {
    fn into_raw(self) -> *mut T;
}

impl<T, E: Display> ResultExt<T, E> for Result<T, E> {
    fn into_raw(self) -> *mut T {
        match self {
            Ok(value) => Box::into_raw(Box::new(value)),
            Err(e) => {
                set_error(e);
                std::ptr::null_mut()
            }
        }
    }
}

pub type Vec2 = glm::TVec2<f32>;
pub type Vec3 = glm::TVec3<f32>;

struct DummyVec2;

impl NormWrapper<DummyVec2> for Vec2 {
    type Norm = f32;

    fn norm(&self) -> Self::Norm {
        self.norm()
    }
}

pub struct DummyVec3;

impl NormWrapper<DummyVec3> for Vec3 {
    type Norm = f32;

    fn norm(&self) -> Self::Norm {
        self.norm()
    }
}

/// A (three-dimensional) ASDF spline.
pub type AsdfPosSpline3 = AsdfPosSpline<f32, Vec3, DummyVec3>;
pub type AsdfCubicCurve3 = PiecewiseCubicCurve<f32, Vec3>;
pub type AsdfCubicCurve2 = PiecewiseCubicCurve<f32, Vec2>;
pub type AsdfCubicCurve1 = PiecewiseCubicCurve<f32, f32>;
pub type AsdfMonotoneCubic = MonotoneCubicSpline<f32>;

/// Creates an `AsdfPosSpline3`.
///
/// Each element in `positions` (3D coordinates) and `tcb`
/// (tension, continuity, bias) contains *three* `float` values,
/// `times` and `speeds` contain one `float` per element.
///
/// # Safety
///
/// All input pointers must be valid for the corresponding `*_count` numbers
/// of elements (not bytes).
#[no_mangle]
pub unsafe extern "C" fn asdf_asdfposspline3(
    positions: *const f32,
    positions_count: size_t,
    times: *const f32,
    times_count: size_t,
    speeds: *const f32,
    speeds_count: size_t,
    tcb: *const f32,
    tcb_count: size_t,
    closed: bool,
) -> *mut AsdfPosSpline3 {
    let positions: Vec<_> = slice::from_raw_parts(positions as *const [f32; 3], positions_count)
        .iter()
        .map(|coords| Vec3::from_column_slice(coords))
        .collect();
    let times: Vec<_> = slice::from_raw_parts(times, times_count)
        .iter()
        .map(|&t| if t.is_nan() { None } else { Some(t) })
        .collect();
    let speeds: Vec<_> = slice::from_raw_parts(speeds, speeds_count)
        .iter()
        .map(|&t| if t.is_nan() { None } else { Some(t) })
        .collect();
    let tcb = slice::from_raw_parts(tcb as *const [f32; 3], tcb_count);
    AsdfPosSpline3::new(&positions, &times, &speeds, tcb, closed).into_raw()
}

/// Frees an `AsdfPosSpline3`
///
/// # Safety
///
/// The pointer must have been obtained with `asdf_asdfposspline3()`.
/// Each pointer can only be freed once.
#[no_mangle]
pub unsafe extern "C" fn asdf_asdfposspline3_free(ptr: *mut AsdfPosSpline3) {
    if !ptr.is_null() {
        Box::from_raw(ptr);
    }
}

/// Returns curve value(s) at given time(s).
///
/// # Safety
///
/// All pointers must be valid.
/// `times` contains one `float` per element,
/// `output` must provide space for *three* `float`s per element.
#[no_mangle]
pub unsafe extern "C" fn asdf_asdfposspline3_evaluate(
    ptr: *mut AsdfPosSpline3,
    times: *const f32,
    count: size_t,
    output: *mut f32,
) {
    assert!(!ptr.is_null());
    let curve = &mut *ptr;
    for i in 0..count {
        let v = curve.evaluate(*times.add(i));
        *output.add(3 * i) = *v.as_ptr();
        *output.add(3 * i + 1) = *v.as_ptr().add(1);
        *output.add(3 * i + 2) = *v.as_ptr().add(2);
    }
}

/// Provides a pointer to (and number of) grid elements.
///
/// # Safety
///
/// All pointers must be valid.
#[no_mangle]
pub unsafe extern "C" fn asdf_asdfposspline3_grid(
    ptr: *mut AsdfPosSpline3,
    output: *mut *const f32,
) -> size_t {
    assert!(!ptr.is_null());
    let curve = &mut *ptr;
    let grid = curve.grid();
    *output = grid.as_ptr();
    grid.len()
}

/// Creates a three-dimensional KB-spline.
///
/// Each element in `positions` (3D coordinates) and `tcb`
/// (tension, continuity, bias) contains *three* `float` values,
///
/// # Safety
///
/// All input pointers must be valid for the corresponding `*_count` numbers
/// of elements (not bytes).
#[no_mangle]
pub unsafe extern "C" fn asdf_centripetalkochanekbartelsspline3(
    positions: *const f32,
    positions_count: size_t,
    tcb: *const f32,
    tcb_count: size_t,
    closed: bool,
) -> *mut AsdfCubicCurve3 {
    let positions = slice::from_raw_parts(positions as *const [f32; 3], positions_count);
    let positions: Vec<_> = positions
        .iter()
        .map(|coords| Vec3::from_column_slice(coords))
        .collect();
    let tcb = slice::from_raw_parts(tcb as *const [f32; 3], tcb_count);
    PiecewiseCubicCurve::new_centripetal_kochanek_bartels(&positions, tcb, closed).into_raw()
}

/// Creates a two-dimensional KB-spline.
///
/// Each element in `positions` (2D coordinates) contains *two* `float` values,
/// each element in `tcb` (tension, continuity, bias) contains *three* `float` values,
///
/// # Safety
///
/// All input pointers must be valid for the corresponding `*_count` numbers
/// of elements (not bytes).
#[no_mangle]
pub unsafe extern "C" fn asdf_centripetalkochanekbartelsspline2(
    positions: *const f32,
    positions_count: size_t,
    tcb: *const f32,
    tcb_count: size_t,
    closed: bool,
) -> *mut AsdfCubicCurve2 {
    let positions = slice::from_raw_parts(positions as *const [f32; 2], positions_count);
    let positions: Vec<_> = positions
        .iter()
        .map(|coords| Vec2::from_column_slice(coords))
        .collect();
    let tcb = slice::from_raw_parts(tcb as *const [f32; 3], tcb_count);
    PiecewiseCubicCurve::new_centripetal_kochanek_bartels(&positions, tcb, closed).into_raw()
}

/// Creates a one-dimensional shape-preserving cubic spline.
///
/// # Safety
///
/// All input pointers must be valid for the corresponding `*_count` numbers
/// of elements (not bytes).
#[no_mangle]
pub unsafe extern "C" fn asdf_shapepreservingcubicspline(
    values: *const f32,
    values_count: size_t,
    grid: *const f32,
    grid_count: size_t,
    closed: bool,
) -> *mut AsdfCubicCurve1 {
    let values = slice::from_raw_parts(values, values_count);
    let grid = slice::from_raw_parts(grid, grid_count);
    PiecewiseCubicCurve::new_shape_preserving(values, grid, closed).into_raw()
}

/// Creates a one-dimensional shape-preserving cubic spline (given values and slopes).
///
/// # Safety
///
/// All input pointers must be valid for the corresponding `*_count` numbers
/// of elements (not bytes).
#[no_mangle]
pub unsafe extern "C" fn asdf_shapepreservingcubicspline_with_slopes(
    values: *const f32,
    values_count: size_t,
    slopes: *const f32,
    slopes_count: size_t,
    grid: *const f32,
    grid_count: size_t,
    closed: bool,
) -> *mut AsdfCubicCurve1 {
    let values = slice::from_raw_parts(values, values_count);
    let slopes = slice::from_raw_parts(slopes, slopes_count);
    let slopes: Vec<_> = slopes
        .iter()
        .map(|&x| if x.is_nan() { None } else { Some(x) })
        .collect();
    let grid = slice::from_raw_parts(grid, grid_count);
    PiecewiseCubicCurve::new_shape_preserving_with_slopes(values, slopes, grid, closed).into_raw()
}

/// Creates a one-dimensional monotone cubic spline.
///
/// # Safety
///
/// All input pointers must be valid for the corresponding `*_count` numbers
/// of elements (not bytes).
#[no_mangle]
pub unsafe extern "C" fn asdf_monotonecubic(
    values: *const f32,
    values_count: size_t,
    grid: *const f32,
    grid_count: size_t,
) -> *mut AsdfMonotoneCubic {
    let values = slice::from_raw_parts(values, values_count);
    let grid = slice::from_raw_parts(grid, grid_count);
    MonotoneCubicSpline::new(values, grid).into_raw()
}

/// Creates a one-dimensional monotone cubic spline (given values and slopes).
///
/// # Safety
///
/// All input pointers must be valid for the corresponding `*_count` numbers
/// of elements (not bytes).
#[no_mangle]
pub unsafe extern "C" fn asdf_monotonecubic_with_slopes(
    values: *const f32,
    values_count: size_t,
    slopes: *const f32,
    slopes_count: size_t,
    grid: *const f32,
    grid_count: size_t,
) -> *mut AsdfMonotoneCubic {
    let values = slice::from_raw_parts(values, values_count);
    let slopes = slice::from_raw_parts(slopes, slopes_count);
    let slopes: Vec<_> = slopes
        .iter()
        .map(|&x| if x.is_nan() { None } else { Some(x) })
        .collect();
    let grid = slice::from_raw_parts(grid, grid_count);
    MonotoneCubicSpline::with_slopes(values, slopes, grid).into_raw()
}

/// Frees an `AsdfMonotoneCubic`
///
/// # Safety
///
/// The pointer must have been obtained with `asdf_monotonecubic()`.
/// Each pointer can only be freed once.
#[no_mangle]
pub unsafe extern "C" fn asdf_monotonecubic_free(ptr: *mut AsdfMonotoneCubic) {
    if !ptr.is_null() {
        Box::from_raw(ptr);
    }
}

/// Returns a pointer to `AsdfCubicCurve1` from `AsdfMonotoneCubic`.
///
/// # Safety
///
/// The pointer must be valid.
#[no_mangle]
pub unsafe extern "C" fn asdf_monotonecubic_inner(
    ptr: *mut AsdfMonotoneCubic,
) -> *const AsdfCubicCurve1 {
    assert!(!ptr.is_null());
    let curve = &mut *ptr;
    curve.inner_ref()
}

/// Returns the time instance(s) for the given value(s).
///
/// # Safety
///
/// All pointers must be valid.
/// `values` contains one `float` per element,
/// `output` must provide space for one `float` per element.
#[no_mangle]
pub unsafe extern "C" fn asdf_monotonecubic_get_time(
    ptr: *mut AsdfMonotoneCubic,
    values: *const f32,
    count: size_t,
    output: *mut f32,
) {
    assert!(!ptr.is_null());
    let curve = &mut *ptr;
    for i in 0..count {
        *output.add(i) = curve.get_time(*values.add(i)).unwrap_or(std::f32::NAN);
    }
}

// TODO: avoid duplication for 1, 2 and 3 dimensions ...

/// Frees an `AsdfCubicCurve3`
///
/// # Safety
///
/// The pointer must have been obtained with `asdf_centripetalkochanekbartelsspline3()`.
/// Each pointer can only be freed once.
#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve3_free(ptr: *mut AsdfCubicCurve3) {
    if !ptr.is_null() {
        Box::from_raw(ptr);
    }
}

/// Returns curve value(s) at given time(s).
///
/// # Safety
///
/// All pointers must be valid.
/// `times` contains one `float` per element,
/// `output` must provide space for *three* `float`s per element.
#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve3_evaluate(
    ptr: *mut AsdfCubicCurve3,
    times: *const f32,
    count: size_t,
    output: *mut f32,
) {
    assert!(!ptr.is_null());
    let curve = &mut *ptr;
    for i in 0..count {
        let v = curve.evaluate(*times.add(i));
        *output.add(3 * i) = *v.as_ptr();
        *output.add(3 * i + 1) = *v.as_ptr().add(1);
        *output.add(3 * i + 2) = *v.as_ptr().add(2);
    }
}

/// Provides a pointer to (and number of) grid elements.
///
/// # Safety
///
/// All pointers must be valid.
#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve3_grid(
    ptr: *mut AsdfCubicCurve3,
    output: *mut *const f32,
) -> size_t {
    assert!(!ptr.is_null());
    let curve = &mut *ptr;
    let grid = curve.grid();
    *output = grid.as_ptr();
    grid.len()
}

/// Frees an `AsdfCubicCurve2`
///
/// # Safety
///
/// The pointer must have been obtained with `asdf_centripetalkochanekbartelsspline2()`.
/// Each pointer can only be freed once.
#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve2_free(ptr: *mut AsdfCubicCurve2) {
    if !ptr.is_null() {
        Box::from_raw(ptr);
    }
}

/// Returns curve value(s) at given time(s).
///
/// # Safety
///
/// All pointers must be valid.
/// `times` contains one `float` per element,
/// `output` must provide space for *two* `float`s per element.
#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve2_evaluate(
    ptr: *mut AsdfCubicCurve2,
    times: *const f32,
    count: size_t,
    output: *mut f32,
) {
    assert!(!ptr.is_null());
    let curve = &mut *ptr;
    for i in 0..count {
        let v = curve.evaluate(*times.add(i));
        *output.add(2 * i) = *v.as_ptr();
        *output.add(2 * i + 1) = *v.as_ptr().add(1);
    }
}

/// Provides a pointer to (and number of) grid elements.
///
/// # Safety
///
/// All pointers must be valid.
#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve2_grid(
    ptr: *mut AsdfCubicCurve2,
    output: *mut *const f32,
) -> size_t {
    assert!(!ptr.is_null());
    let curve = &mut *ptr;
    let grid = curve.grid();
    *output = grid.as_ptr();
    grid.len()
}

/// Frees an `AsdfCubicCurve1`
///
/// # Safety
///
/// The pointer must have been obtained with `asdf_shapepreservingcubicspline()` or
/// `asdf_shapepreservingcubicspline_with_slopes()`.
/// Each pointer can only be freed once.
#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve1_free(ptr: *mut AsdfCubicCurve1) {
    if !ptr.is_null() {
        Box::from_raw(ptr);
    }
}

/// Returns curve value(s) at given time(s).
///
/// # Safety
///
/// All pointers must be valid.
/// `times` contains one `float` per element,
/// `output` must provide space for *one* `float` per element.
#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve1_evaluate(
    ptr: *mut AsdfCubicCurve1,
    times: *const f32,
    count: size_t,
    output: *mut f32,
) {
    assert!(!ptr.is_null());
    let curve = &mut *ptr;
    for i in 0..count {
        *output.add(i) = curve.evaluate(*times.add(i));
    }
}

/// Provides a pointer to (and number of) grid elements.
///
/// # Safety
///
/// All pointers must be valid.
#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve1_grid(
    ptr: *mut AsdfCubicCurve1,
    output: *mut *const f32,
) -> size_t {
    assert!(!ptr.is_null());
    let curve = &mut *ptr;
    let grid = curve.grid();
    *output = grid.as_ptr();
    grid.len()
}
