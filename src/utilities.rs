use num_traits::{one, zero, Float, FromPrimitive};

/// <https://en.wikipedia.org/wiki/Bisection_method>
///
/// Root must be within `[xmin, xmax]`, otherwise one of those is returned
/// (whichever has a function value closer to zero).
pub fn bisect<T, F>(f: F, mut xmin: T, mut xmax: T, xtol: T, max_calls: usize) -> T
where
    T: Float,
    F: Fn(T) -> T,
{
    assert!(xmin <= xmax);
    let two = one::<T>() + one();
    let mut calls: usize = 0;
    let mut fmin = f(xmin);
    calls += 1;
    if fmin == zero() {
        return xmin;
    }
    let mut fmax = f(xmax);
    calls += 1;
    if fmax == zero() {
        return xmax;
    }
    assert!(max_calls >= calls);
    if fmin * fmax < zero() {
        while (max_calls - calls) > 0 && (xmax - xmin) > xtol {
            let xmid = (xmin + xmax) / two;
            if xmid == xmin || xmid == xmax {
                return xmid;
            }
            let fmid = f(xmid);
            calls += 1;
            if fmid == zero() {
                return xmid;
            }
            if fmin * fmid < zero() {
                xmax = xmid;
                fmax = fmid;
            } else {
                xmin = xmid;
                fmin = fmid;
            }
        }
    }
    if fmin.abs() < fmax.abs() {
        xmin
    } else {
        xmax
    }
    // TODO: return number of calls?
    // TODO: return function value that's supposedly zero?
}

/// Gauss-Legendre quadrature of order 13.
///
/// https://en.wikipedia.org/wiki/Gaussian_quadrature
///
/// Arrays were generated with scipy.special.roots_legendre(13).
/// 13th order typically leads to results within single-precision
/// accuracy [citation needed].
///
/// If used with `f64`, this will still only have `f32` accuracy!
///
/// See also https://pomax.github.io/bezierinfo/legendre-gauss.html
pub fn gauss_legendre13<T: Float, F>(f: F, a: T, b: T) -> T
where
    T: Float + FromPrimitive,
    F: Fn(T) -> T,
{
    // TODO: separate versions for f32 and f64?
    #[allow(clippy::unreadable_literal)]
    let times = [
        T::from_f64(-0.9841830547185881).unwrap(),
        T::from_f64(-0.9175983992229779).unwrap(),
        T::from_f64(-0.8015780907333099).unwrap(),
        T::from_f64(-0.6423493394403403).unwrap(),
        T::from_f64(-0.44849275103644687).unwrap(),
        T::from_f64(-0.23045831595513483).unwrap(),
        T::from_f64(0.0).unwrap(),
        T::from_f64(0.23045831595513483).unwrap(),
        T::from_f64(0.44849275103644687).unwrap(),
        T::from_f64(0.6423493394403403).unwrap(),
        T::from_f64(0.8015780907333099).unwrap(),
        T::from_f64(0.9175983992229779).unwrap(),
        T::from_f64(0.9841830547185881).unwrap(),
    ];
    #[allow(clippy::unreadable_literal)]
    let weights = [
        T::from_f64(0.04048400476531615).unwrap(),
        T::from_f64(0.0921214998377276).unwrap(),
        T::from_f64(0.1388735102197876).unwrap(),
        T::from_f64(0.17814598076194554).unwrap(),
        T::from_f64(0.20781604753688862).unwrap(),
        T::from_f64(0.2262831802628975).unwrap(),
        T::from_f64(0.23255155323087406).unwrap(),
        T::from_f64(0.2262831802628975).unwrap(),
        T::from_f64(0.20781604753688862).unwrap(),
        T::from_f64(0.17814598076194554).unwrap(),
        T::from_f64(0.1388735102197876).unwrap(),
        T::from_f64(0.0921214998377276).unwrap(),
        T::from_f64(0.04048400476531615).unwrap(),
    ];
    assert_eq!(times.len(), weights.len());
    let two = one::<T>() + one();
    let sum = (0..times.len())
        .map(|i| weights[i] * f((b - a) * times[i] / two + (a + b) / two))
        .fold(zero(), |acc, x| acc + x);
    (b - a) * sum / two
}

pub enum GridError {
    GridNan { index: usize },
    GridNotAscending { index: usize },
}

pub fn check_grid<S: crate::Scalar>(grid: &[S]) -> Result<(), GridError> {
    use GridError::*;
    if let Some(index) = grid.iter().copied().position(S::is_nan) {
        return Err(GridNan { index });
    }
    if let Some(index) = grid.windows(2).position(|w| w[0] >= w[1]) {
        return Err(GridNotAscending { index: index + 1 });
    }
    Ok(())
}
