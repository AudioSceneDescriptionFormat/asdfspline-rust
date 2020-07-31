/// <https://en.wikipedia.org/wiki/Bisection_method>
///
/// Root must be within `[xmin, xmax]`, otherwise one of those is returned
/// (whichever has a function value closer to zero).
pub fn bisect<F>(f: F, mut xmin: f32, mut xmax: f32, xtol: f32, max_calls: usize) -> f32
where
    F: Fn(f32) -> f32,
{
    assert!(xmin <= xmax);
    let mut calls: usize = 0;
    let mut fmin = f(xmin);
    calls += 1;
    if fmin == 0.0 {
        return xmin;
    }
    let mut fmax = f(xmax);
    calls += 1;
    if fmax == 0.0 {
        return xmax;
    }
    assert!(max_calls >= calls);
    if fmin * fmax < 0.0 {
        while (max_calls - calls) > 0 && (xmax - xmin) > xtol {
            let xmid = (xmin + xmax) / 2.0;
            if xmid <= xmin || xmid >= xmax {
                return xmid;
            }
            let fmid = f(xmid);
            calls += 1;
            if fmid == 0.0 {
                return xmid;
            }
            if fmin * fmid < 0.0 {
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
/// See also https://pomax.github.io/bezierinfo/legendre-gauss.html
pub fn gauss_legendre13<F>(f: F, a: f32, b: f32) -> f32
where
    F: Fn(f32) -> f32,
{
    #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
    let times = [
        -0.9841830547185881,
        -0.9175983992229779,
        -0.8015780907333099,
        -0.6423493394403403,
        -0.44849275103644687,
        -0.23045831595513483,
        0.0,
        0.23045831595513483,
        0.44849275103644687,
        0.6423493394403403,
        0.8015780907333099,
        0.9175983992229779,
        0.9841830547185881,
    ];
    #[allow(clippy::unreadable_literal, clippy::excessive_precision)]
    let weights = [
        0.04048400476531615,
        0.0921214998377276,
        0.1388735102197876,
        0.17814598076194554,
        0.20781604753688862,
        0.2262831802628975,
        0.23255155323087406,
        0.2262831802628975,
        0.20781604753688862,
        0.17814598076194554,
        0.1388735102197876,
        0.0921214998377276,
        0.04048400476531615,
    ];
    assert_eq!(times.len(), weights.len());
    let sum = (0..times.len())
        .map(|i| weights[i] * f((b - a) * times[i] / 2.0 + (a + b) / 2.0))
        .fold(0.0, |acc, x| acc + x);
    (b - a) * sum / 2.0
}

pub enum GridError {
    GridNan { index: usize },
    GridNotAscending { index: usize },
}

pub fn check_grid(grid: &[f32]) -> Result<(), GridError> {
    use GridError::*;
    if let Some(index) = grid.iter().copied().position(f32::is_nan) {
        return Err(GridNan { index });
    }
    if let Some(index) = grid.windows(2).position(|w| w[0] >= w[1]) {
        return Err(GridNotAscending { index: index + 1 });
    }
    Ok(())
}
