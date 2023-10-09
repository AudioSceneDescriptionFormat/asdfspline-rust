#![deny(unsafe_op_in_unsafe_fn)]

use std::cell::RefCell;
use std::ffi::CString;
use std::fmt::Display;
use std::mem::MaybeUninit;
use std::slice;

use libc::{c_char, size_t};
use nalgebra::{Vector2, Vector3};

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

trait ResultExt<T, E> {
    fn into_box(self) -> Option<Box<T>>;
}

impl<T, E: Display> ResultExt<T, E> for Result<T, E> {
    fn into_box(self) -> Option<Box<T>> {
        self.map(Box::new).map_err(|e| set_error(e)).ok()
    }
}

pub type Vec2 = Vector2<f32>;
pub type Vec3 = Vector3<f32>;

pub struct Norm3;

impl NormWrapper<Norm3> for Vec3 {
    fn norm(&self) -> f32 {
        self.norm()
    }
}

/// A (three-dimensional) ASDF spline.
pub type AsdfPosSpline3 = AsdfPosSpline<Vec3, Norm3>;
pub type AsdfCubicCurve3 = PiecewiseCubicCurve<Vec3>;
pub type AsdfCubicCurve2 = PiecewiseCubicCurve<Vec2>;
pub type AsdfCubicCurve1 = PiecewiseCubicCurve<f32>;
pub type AsdfMonotoneCubic = MonotoneCubicSpline;

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
) -> Option<Box<AsdfPosSpline3>> {
    let positions: Vec<_> =
        unsafe { slice::from_raw_parts(positions.cast::<[f32; 3]>(), positions_count) }
            .iter()
            .map(|coords| Vec3::from_column_slice(coords))
            .collect();
    let times: Vec<_> = unsafe { slice::from_raw_parts(times, times_count) }
        .iter()
        .map(|&t| if t.is_nan() { None } else { Some(t) })
        .collect();
    let speeds: Vec<_> = unsafe { slice::from_raw_parts(speeds, speeds_count) }
        .iter()
        .map(|&t| if t.is_nan() { None } else { Some(t) })
        .collect();
    let tcb = unsafe { slice::from_raw_parts(tcb.cast::<[f32; 3]>(), tcb_count) };
    AsdfPosSpline3::new(positions, times, speeds, tcb, closed).into_box()
}

/// Frees an `AsdfPosSpline3`
///
/// # Safety
///
/// The pointer must have been obtained with `asdf_asdfposspline3()`.
/// Each pointer can only be freed once.
/// Passing NULL is allowed.
#[no_mangle]
pub unsafe extern "C" fn asdf_asdfposspline3_free(_: Option<Box<AsdfPosSpline3>>) {}

/// Returns curve value(s) at given time(s).
///
/// # Safety
///
/// All pointers must be valid.
/// `times` contains one `float` per element,
/// `output` must provide space for *three* `float`s per element.
#[no_mangle]
pub unsafe extern "C" fn asdf_asdfposspline3_evaluate(
    curve: &mut AsdfPosSpline3,
    times: *const f32,
    count: size_t,
    output: *mut f32,
) {
    let times = unsafe { std::slice::from_raw_parts(times, count) };
    let output =
        unsafe { std::slice::from_raw_parts_mut(output.cast::<MaybeUninit<[f32; 3]>>(), count) };
    for (time, out) in times.iter().zip(output) {
        *out = MaybeUninit::new(curve.evaluate(*time).into());
    }
}

/// Provides a pointer to (and number of) grid elements.
///
/// # Safety
///
/// All pointers must be valid.
#[no_mangle]
pub unsafe extern "C" fn asdf_asdfposspline3_grid(
    curve: &mut AsdfPosSpline3,
    output: *mut *const f32,
) -> size_t {
    let grid = curve.grid();
    unsafe { output.write(grid.as_ptr()) };
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
) -> Option<Box<AsdfCubicCurve3>> {
    let positions = unsafe { slice::from_raw_parts(positions.cast::<[f32; 3]>(), positions_count) };
    let positions: Vec<_> = positions
        .iter()
        .map(|coords| Vec3::from_column_slice(coords))
        .collect();
    let tcb = unsafe { slice::from_raw_parts(tcb.cast::<[f32; 3]>(), tcb_count) };
    PiecewiseCubicCurve::new_centripetal_kochanek_bartels(&positions, tcb, closed, Vec3::norm)
        .into_box()
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
) -> Option<Box<AsdfCubicCurve2>> {
    let positions = unsafe { slice::from_raw_parts(positions.cast::<[f32; 2]>(), positions_count) };
    let positions: Vec<_> = positions
        .iter()
        .map(|coords| Vec2::from_column_slice(coords))
        .collect();
    let tcb = unsafe { slice::from_raw_parts(tcb.cast::<[f32; 3]>(), tcb_count) };
    PiecewiseCubicCurve::new_centripetal_kochanek_bartels(&positions, tcb, closed, Vec2::norm)
        .into_box()
}

/// Creates a one-dimensional piecewise monotone cubic spline.
///
/// # Safety
///
/// All input pointers must be valid for the corresponding `*_count` numbers
/// of elements (not bytes).
#[no_mangle]
pub unsafe extern "C" fn asdf_piecewisemonotonecubicspline(
    values: *const f32,
    values_count: size_t,
    grid: *const f32,
    grid_count: size_t,
    closed: bool,
) -> Option<Box<AsdfCubicCurve1>> {
    let values = unsafe { slice::from_raw_parts(values, values_count) };
    let grid = unsafe { slice::from_raw_parts(grid, grid_count) };
    PiecewiseCubicCurve::new_piecewise_monotone(values, grid, closed).into_box()
}

/// Creates a one-dimensional piecewise monotone cubic spline (given values and slopes).
///
/// # Safety
///
/// All input pointers must be valid for the corresponding `*_count` numbers
/// of elements (not bytes).
#[no_mangle]
pub unsafe extern "C" fn asdf_piecewisemonotonecubicspline_with_slopes(
    values: *const f32,
    values_count: size_t,
    slopes: *const f32,
    slopes_count: size_t,
    grid: *const f32,
    grid_count: size_t,
    closed: bool,
) -> Option<Box<AsdfCubicCurve1>> {
    let values = unsafe { slice::from_raw_parts(values, values_count) };
    let slopes = unsafe { slice::from_raw_parts(slopes, slopes_count) };
    let slopes: Vec<_> = slopes
        .iter()
        .map(|&x| if x.is_nan() { None } else { Some(x) })
        .collect();
    let grid = unsafe { slice::from_raw_parts(grid, grid_count) };
    PiecewiseCubicCurve::new_piecewise_monotone_with_slopes(values, slopes, grid, closed).into_box()
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
    cyclic: bool,
) -> Option<Box<AsdfMonotoneCubic>> {
    let values = unsafe { slice::from_raw_parts(values, values_count) };
    let grid = unsafe { slice::from_raw_parts(grid, grid_count) };
    MonotoneCubicSpline::new(values, grid, cyclic).into_box()
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
    cyclic: bool,
) -> Option<Box<AsdfMonotoneCubic>> {
    let values = unsafe { slice::from_raw_parts(values, values_count) };
    let slopes = unsafe { slice::from_raw_parts(slopes, slopes_count) };
    let slopes: Vec<_> = slopes
        .iter()
        .map(|&x| if x.is_nan() { None } else { Some(x) })
        .collect();
    let grid = unsafe { slice::from_raw_parts(grid, grid_count) };
    MonotoneCubicSpline::with_slopes(values, slopes, grid, cyclic).into_box()
}

/// Frees an `AsdfMonotoneCubic`
///
/// # Safety
///
/// The pointer must have been obtained with `asdf_monotonecubic()`.
/// Each pointer can only be freed once.
/// Passing NULL is allowed.
#[no_mangle]
pub unsafe extern "C" fn asdf_monotonecubic_free(_: Option<Box<AsdfMonotoneCubic>>) {}

/// Returns a pointer to `AsdfCubicCurve1` from `AsdfMonotoneCubic`.
///
/// # Safety
///
/// The pointer must be valid.
#[no_mangle]
pub unsafe extern "C" fn asdf_monotonecubic_inner(
    curve: &mut AsdfMonotoneCubic,
) -> *const AsdfCubicCurve1 {
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
    curve: &mut AsdfMonotoneCubic,
    values: *const f32,
    count: size_t,
    output: *mut f32,
) {
    let values = unsafe { slice::from_raw_parts(values, count) };
    let output = unsafe { slice::from_raw_parts_mut(output.cast::<MaybeUninit<_>>(), count) };
    for (val, out) in values.iter().zip(output) {
        *out = MaybeUninit::new(curve.get_time(*val).unwrap_or(std::f32::NAN));
    }
}

// TODO: avoid duplication for 1, 2 and 3 dimensions ...

/// Frees an `AsdfCubicCurve3`
///
/// # Safety
///
/// The pointer must have been obtained with `asdf_centripetalkochanekbartelsspline3()`.
/// Each pointer can only be freed once.
/// Passing NULL is allowed.
#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve3_free(_: Option<Box<AsdfCubicCurve3>>) {}

/// Returns curve value(s) at given time(s).
///
/// # Safety
///
/// All pointers must be valid.
/// `times` contains one `float` per element,
/// `output` must provide space for *three* `float`s per element.
#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve3_evaluate(
    curve: &mut AsdfCubicCurve3,
    times: *const f32,
    count: size_t,
    output: *mut f32,
) {
    let times = unsafe { std::slice::from_raw_parts(times, count) };
    let output =
        unsafe { std::slice::from_raw_parts_mut(output.cast::<MaybeUninit<[f32; 3]>>(), count) };
    for (time, out) in times.iter().zip(output) {
        *out = MaybeUninit::new(curve.evaluate(*time).into());
    }
}

/// Provides a pointer to (and number of) grid elements.
///
/// # Safety
///
/// All pointers must be valid.
#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve3_grid(
    curve: &mut AsdfCubicCurve3,
    output: *mut *const f32,
) -> size_t {
    let grid = curve.grid();
    unsafe {
        output.write(grid.as_ptr());
    }
    grid.len()
}

/// Frees an `AsdfCubicCurve2`
///
/// # Safety
///
/// The pointer must have been obtained with `asdf_centripetalkochanekbartelsspline2()`.
/// Each pointer can only be freed once.
#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve2_free(_: Option<Box<AsdfCubicCurve2>>) {}

/// Returns curve value(s) at given time(s).
///
/// # Safety
///
/// All pointers must be valid.
/// `times` contains one `float` per element,
/// `output` must provide space for *two* `float`s per element.
#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve2_evaluate(
    curve: &mut AsdfCubicCurve2,
    times: *const f32,
    count: size_t,
    output: *mut f32,
) {
    let times = unsafe { std::slice::from_raw_parts(times, count) };
    let output =
        unsafe { std::slice::from_raw_parts_mut(output.cast::<MaybeUninit<[f32; 2]>>(), count) };
    for (time, out) in times.iter().zip(output) {
        *out = MaybeUninit::new(curve.evaluate(*time).into());
    }
}

/// Provides a pointer to (and number of) grid elements.
///
/// # Safety
///
/// All pointers must be valid.
#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve2_grid(
    curve: &mut AsdfCubicCurve2,
    output: *mut *const f32,
) -> size_t {
    let grid = curve.grid();
    unsafe {
        output.write(grid.as_ptr());
    }
    grid.len()
}

/// Frees an `AsdfCubicCurve1`
///
/// # Safety
///
/// The pointer must have been obtained with `asdf_piecewisemonotonecubicspline()` or
/// `asdf_piecewisemonotonecubicspline_with_slopes()`.
/// Each pointer can only be freed once.
/// Passing NULL is allowed.
#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve1_free(_: Option<Box<AsdfCubicCurve1>>) {}

/// Returns curve value(s) at given time(s).
///
/// # Safety
///
/// All pointers must be valid.
/// `times` contains one `float` per element,
/// `output` must provide space for *one* `float` per element.
#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve1_evaluate(
    curve: &mut AsdfCubicCurve1,
    times: *const f32,
    count: size_t,
    output: *mut f32,
) {
    let times = unsafe { std::slice::from_raw_parts(times, count) };
    let output =
        unsafe { std::slice::from_raw_parts_mut(output.cast::<MaybeUninit<f32>>(), count) };
    for (time, out) in times.iter().zip(output) {
        *out = MaybeUninit::new(curve.evaluate(*time));
    }
}

/// Provides a pointer to (and number of) grid elements.
///
/// # Safety
///
/// All pointers must be valid.
#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve1_grid(
    curve: &mut AsdfCubicCurve1,
    output: *mut *const f32,
) -> size_t {
    let grid = curve.grid();
    unsafe {
        output.write(grid.as_ptr());
    }
    grid.len()
}
