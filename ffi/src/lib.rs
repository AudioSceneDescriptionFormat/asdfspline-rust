use std::cell::RefCell;
use std::ffi::CString;
use std::fmt::Display;
use std::panic::{catch_unwind, UnwindSafe};
use std::slice;

use libc::{c_char, size_t};
use nalgebra_glm as glm;

use asdfspline::{
    make_centripetal_kochanek_bartels_spline, make_shape_preserving_cubic_spline,
    make_shape_preserving_cubic_spline_with_slopes, MonotoneCubicSpline, PiecewiseCubicCurve,
};

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

fn handle_errors<F, T>(f: F, optb: T) -> T
where
    F: FnOnce() -> T + UnwindSafe,
{
    match catch_unwind(f) {
        Ok(value) => value,
        Err(e) => {
            if let Some(e) = e.downcast_ref::<&str>() {
                set_error(*e);
            } else if let Some(e) = e.downcast_ref::<String>() {
                set_error(e);
            } else {
                set_error("unknown error");
            }
            optb
        }
    }
}

trait ResultExt<T, E: Display> {
    fn unwrap_display(self) -> T;
}

impl<T, E: Display> ResultExt<T, E> for Result<T, E> {
    fn unwrap_display(self) -> T {
        match self {
            Ok(value) => value,
            Err(e) => panic!(e.to_string()),
        }
    }
}

pub type Vec2 = glm::TVec2<f32>;
pub type Vec3 = glm::TVec3<f32>;
pub type AsdfSpline = asdfspline::AsdfSpline<f32, Vec3>;
pub type AsdfCubicCurve3 = PiecewiseCubicCurve<f32, Vec3>;
pub type AsdfCubicCurve2 = PiecewiseCubicCurve<f32, Vec2>;
pub type AsdfCubicCurve1 = PiecewiseCubicCurve<f32, f32>;
pub type AsdfMonotoneCubic = MonotoneCubicSpline<f32>;

#[no_mangle]
pub unsafe extern "C" fn asdf_asdfspline(
    positions: *const f32,
    positions_count: size_t,
    times: *const f32,
    times_count: size_t,
    speeds: *const f32,
    speeds_count: size_t,
    tcb: *const f32,
    tcb_count: size_t,
    closed: bool,
) -> *mut AsdfSpline {
    handle_errors(
        || {
            let positions: Vec<_> =
                slice::from_raw_parts(positions as *const [f32; 3], positions_count)
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
            let curve = AsdfSpline::new(positions, times, speeds, tcb, closed, |v| v.norm())
                .unwrap_display();
            Box::into_raw(Box::new(curve))
        },
        std::ptr::null_mut(),
    )
}

#[no_mangle]
pub unsafe extern "C" fn asdf_asdfspline_free(ptr: *mut AsdfSpline) {
    if !ptr.is_null() {
        Box::from_raw(ptr);
    }
}

#[no_mangle]
pub unsafe extern "C" fn asdf_asdfspline_evaluate(
    ptr: *mut AsdfSpline,
    times: *const f32,
    count: size_t,
    output: *mut f32,
) {
    handle_errors(
        || {
            assert!(!ptr.is_null());
            let curve = &mut *ptr;
            for i in 0..count {
                let v = curve.evaluate(*times.add(i), |v| v.norm());
                *output.add(3 * i) = *v.as_ptr();
                *output.add(3 * i + 1) = *v.as_ptr().add(1);
                *output.add(3 * i + 2) = *v.as_ptr().add(2);
            }
        },
        (),
    )
}

#[no_mangle]
pub unsafe extern "C" fn asdf_asdfspline_grid(
    ptr: *mut AsdfSpline,
    output: *mut *const f32,
) -> size_t {
    handle_errors(
        || {
            assert!(!ptr.is_null());
            let curve = &mut *ptr;
            let grid = curve.grid();
            *output = grid.as_ptr();
            grid.len()
        },
        0,
    )
}

#[no_mangle]
pub unsafe extern "C" fn asdf_centripetalkochanekbartelsspline3(
    vertices: *const f32,
    vertices_count: size_t,
    tcb: *const f32,
    tcb_count: size_t,
    closed: bool,
) -> *mut AsdfCubicCurve3 {
    handle_errors(
        || {
            let vertices = slice::from_raw_parts(vertices as *const [f32; 3], vertices_count);
            let vertices: Vec<_> = vertices
                .iter()
                .map(|coords| Vec3::from_column_slice(coords))
                .collect();
            let tcb = slice::from_raw_parts(tcb as *const [f32; 3], tcb_count);
            let curve =
                make_centripetal_kochanek_bartels_spline(vertices, tcb, closed, |v| v.norm())
                    .unwrap_display();
            Box::into_raw(Box::new(curve))
        },
        std::ptr::null_mut(),
    )
}

#[no_mangle]
pub unsafe extern "C" fn asdf_centripetalkochanekbartelsspline2(
    vertices: *const f32,
    vertices_count: size_t,
    tcb: *const f32,
    tcb_count: size_t,
    closed: bool,
) -> *mut AsdfCubicCurve2 {
    handle_errors(
        || {
            let vertices = slice::from_raw_parts(vertices as *const [f32; 2], vertices_count);
            let vertices: Vec<_> = vertices
                .iter()
                .map(|coords| Vec2::from_column_slice(coords))
                .collect();
            let tcb = slice::from_raw_parts(tcb as *const [f32; 3], tcb_count);
            let curve =
                make_centripetal_kochanek_bartels_spline(vertices, tcb, closed, |v| v.norm())
                    .unwrap_display();
            Box::into_raw(Box::new(curve))
        },
        std::ptr::null_mut(),
    )
}

#[no_mangle]
pub unsafe extern "C" fn asdf_shapepreservingcubicspline(
    values: *const f32,
    values_count: size_t,
    grid: *const f32,
    grid_count: size_t,
    closed: bool,
) -> *mut AsdfCubicCurve1 {
    handle_errors(
        || {
            let values = slice::from_raw_parts(values, values_count);
            let grid = slice::from_raw_parts(grid, grid_count);
            let curve = make_shape_preserving_cubic_spline(values, grid, closed).unwrap_display();
            Box::into_raw(Box::new(curve))
        },
        std::ptr::null_mut(),
    )
}

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
    handle_errors(
        || {
            let values = slice::from_raw_parts(values, values_count);
            let slopes = slice::from_raw_parts(slopes, slopes_count);
            let slopes: Vec<_> = slopes
                .iter()
                .map(|&x| if x.is_nan() { None } else { Some(x) })
                .collect();
            let grid = slice::from_raw_parts(grid, grid_count);
            let curve =
                make_shape_preserving_cubic_spline_with_slopes(values, slopes, grid, closed)
                    .unwrap_display();
            Box::into_raw(Box::new(curve))
        },
        std::ptr::null_mut(),
    )
}

#[no_mangle]
pub unsafe extern "C" fn asdf_monotonecubic(
    values: *const f32,
    values_count: size_t,
    grid: *const f32,
    grid_count: size_t,
) -> *mut AsdfMonotoneCubic {
    handle_errors(
        || {
            let values = slice::from_raw_parts(values, values_count);
            let grid = slice::from_raw_parts(grid, grid_count);
            let curve = MonotoneCubicSpline::new(values, grid).unwrap_display();
            Box::into_raw(Box::new(curve))
        },
        std::ptr::null_mut(),
    )
}

#[no_mangle]
pub unsafe extern "C" fn asdf_monotonecubic_with_slopes(
    values: *const f32,
    values_count: size_t,
    slopes: *const f32,
    slopes_count: size_t,
    grid: *const f32,
    grid_count: size_t,
) -> *mut AsdfMonotoneCubic {
    handle_errors(
        || {
            let values = slice::from_raw_parts(values, values_count);
            let slopes = slice::from_raw_parts(slopes, slopes_count);
            let slopes: Vec<_> = slopes
                .iter()
                .map(|&x| if x.is_nan() { None } else { Some(x) })
                .collect();
            let grid = slice::from_raw_parts(grid, grid_count);
            let curve = MonotoneCubicSpline::with_slopes(values, slopes, grid).unwrap_display();
            Box::into_raw(Box::new(curve))
        },
        std::ptr::null_mut(),
    )
}

#[no_mangle]
pub unsafe extern "C" fn asdf_monotonecubic_free(ptr: *mut AsdfMonotoneCubic) {
    if !ptr.is_null() {
        Box::from_raw(ptr);
    }
}

#[no_mangle]
pub unsafe extern "C" fn asdf_monotonecubic_inner(
    ptr: *mut AsdfMonotoneCubic,
) -> *const AsdfCubicCurve1 {
    handle_errors(
        || {
            assert!(!ptr.is_null());
            let curve = &mut *ptr;
            curve.inner_ref()
        },
        std::ptr::null(),
    )
}

#[no_mangle]
pub unsafe extern "C" fn asdf_monotonecubic_get_time(
    ptr: *mut AsdfMonotoneCubic,
    values: *const f32,
    count: size_t,
    output: *mut f32,
) {
    handle_errors(
        || {
            assert!(!ptr.is_null());
            let curve = &mut *ptr;
            for i in 0..count {
                *output.add(i) = curve.get_time(*values.add(i)).unwrap_or(std::f32::NAN);
            }
        },
        (),
    )
}

// TODO: avoid duplication for 1, 2 and 3 dimensions ...

#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve3_free(ptr: *mut AsdfCubicCurve3) {
    if !ptr.is_null() {
        Box::from_raw(ptr);
    }
}

#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve3_evaluate(
    ptr: *mut AsdfCubicCurve3,
    times: *const f32,
    count: size_t,
    output: *mut f32,
) {
    handle_errors(
        || {
            assert!(!ptr.is_null());
            let curve = &mut *ptr;
            for i in 0..count {
                let v = curve.evaluate(*times.add(i));
                *output.add(3 * i) = *v.as_ptr();
                *output.add(3 * i + 1) = *v.as_ptr().add(1);
                *output.add(3 * i + 2) = *v.as_ptr().add(2);
            }
        },
        (),
    )
}

#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve3_grid(
    ptr: *mut AsdfCubicCurve3,
    output: *mut *const f32,
) -> size_t {
    handle_errors(
        || {
            assert!(!ptr.is_null());
            let curve = &mut *ptr;
            let grid = curve.grid();
            *output = grid.as_ptr();
            grid.len()
        },
        0,
    )
}

#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve2_free(ptr: *mut AsdfCubicCurve2) {
    if !ptr.is_null() {
        Box::from_raw(ptr);
    }
}

#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve2_evaluate(
    ptr: *mut AsdfCubicCurve2,
    times: *const f32,
    count: size_t,
    output: *mut f32,
) {
    handle_errors(
        || {
            assert!(!ptr.is_null());
            let curve = &mut *ptr;
            for i in 0..count {
                let v = curve.evaluate(*times.add(i));
                *output.add(2 * i) = *v.as_ptr();
                *output.add(2 * i + 1) = *v.as_ptr().add(1);
            }
        },
        (),
    )
}

#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve2_grid(
    ptr: *mut AsdfCubicCurve2,
    output: *mut *const f32,
) -> size_t {
    handle_errors(
        || {
            assert!(!ptr.is_null());
            let curve = &mut *ptr;
            let grid = curve.grid();
            *output = grid.as_ptr();
            grid.len()
        },
        0,
    )
}

#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve1_free(ptr: *mut AsdfCubicCurve1) {
    if !ptr.is_null() {
        Box::from_raw(ptr);
    }
}

#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve1_evaluate(
    ptr: *mut AsdfCubicCurve1,
    times: *const f32,
    count: size_t,
    output: *mut f32,
) {
    handle_errors(
        || {
            assert!(!ptr.is_null());
            let curve = &mut *ptr;
            for i in 0..count {
                *output.add(i) = curve.evaluate(*times.add(i));
            }
        },
        (),
    )
}

#[no_mangle]
pub unsafe extern "C" fn asdf_cubiccurve1_grid(
    ptr: *mut AsdfCubicCurve1,
    output: *mut *const f32,
) -> size_t {
    handle_errors(
        || {
            assert!(!ptr.is_null());
            let curve = &mut *ptr;
            let grid = curve.grid();
            *output = grid.as_ptr();
            grid.len()
        },
        0,
    )
}
