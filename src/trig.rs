//! Trigonometric functions implemented with Taylor series and range reduction.

use crate::arithmetic::{abs, modulo};
use crate::constants;
use crate::powers::sqrt;

/// Reduce angle to [-π, π].
fn reduce_angle(x: f64) -> f64 {
    let mut r = modulo(x + constants::PI, constants::TAU) - constants::PI;
    if r <= -constants::PI { r += constants::TAU; }
    r
}

/// Sine using Taylor series with range reduction.
pub fn sin(x: f64) -> f64 {
    if x != x { return f64::NAN; }
    let x = reduce_angle(x);
    let x2 = x * x;
    let mut term = x;
    let mut sum = x;
    for i in 1..20 {
        term *= -x2 / ((2 * i) as f64 * (2 * i + 1) as f64);
        sum += term;
        if abs(term) < 1e-16 * abs(sum) { break; }
    }
    sum
}

/// Cosine using Taylor series with range reduction.
pub fn cos(x: f64) -> f64 {
    if x != x { return f64::NAN; }
    let x = reduce_angle(x);
    let x2 = x * x;
    let mut term = 1.0;
    let mut sum = 1.0;
    for i in 1..20 {
        term *= -x2 / ((2 * i - 1) as f64 * (2 * i) as f64);
        sum += term;
        if abs(term) < 1e-16 * abs(sum) { break; }
    }
    sum
}

/// Tangent: sin(x) / cos(x).
pub fn tan(x: f64) -> f64 {
    let c = cos(x);
    if abs(c) < 1e-15 { return if sin(x) > 0.0 { f64::INFINITY } else { f64::NEG_INFINITY }; }
    sin(x) / c
}

/// Cotangent: cos(x) / sin(x).
pub fn cot(x: f64) -> f64 {
    let s = sin(x);
    if abs(s) < 1e-15 { return if cos(x) > 0.0 { f64::INFINITY } else { f64::NEG_INFINITY }; }
    cos(x) / s
}

/// Secant: 1 / cos(x).
pub fn sec(x: f64) -> f64 { 1.0 / cos(x) }

/// Cosecant: 1 / sin(x).
pub fn csc(x: f64) -> f64 { 1.0 / sin(x) }

/// Compute sin and cos simultaneously.
pub fn sincos(x: f64) -> (f64, f64) { (sin(x), cos(x)) }

/// Arcsine using Taylor series + Newton refinement.
pub fn asin(x: f64) -> f64 {
    if x < -1.0 || x > 1.0 { return f64::NAN; }
    if x == 1.0 { return constants::FRAC_PI_2; }
    if x == -1.0 { return -constants::FRAC_PI_2; }
    if abs(x) > 0.7 {
        let sign = if x < 0.0 { -1.0 } else { 1.0 };
        return sign * (constants::FRAC_PI_2 - asin(sqrt(1.0 - x * x)));
    }
    // Taylor series: asin(x) = x + x^3/6 + 3x^5/40 + ...
    let x2 = x * x;
    let mut term = x;
    let mut sum = x;
    for n in 1..40 {
        term *= x2 * (2 * n - 1) as f64 * (2 * n - 1) as f64 / ((2 * n) as f64 * (2 * n + 1) as f64);
        sum += term;
        if abs(term) < 1e-16 { break; }
    }
    sum
}

/// Arccosine: π/2 - asin(x).
pub fn acos(x: f64) -> f64 {
    if x < -1.0 || x > 1.0 { return f64::NAN; }
    constants::FRAC_PI_2 - asin(x)
}

/// Arctangent using range reduction and Taylor series.
pub fn atan(x: f64) -> f64 {
    if x != x { return f64::NAN; }
    if x == f64::INFINITY { return constants::FRAC_PI_2; }
    if x == f64::NEG_INFINITY { return -constants::FRAC_PI_2; }
    // Use symmetry for negative
    if x < 0.0 { return -atan(-x); }
    // Range reduction: for x > 1, use atan(x) = pi/2 - atan(1/x)
    if x > 1.0 { return constants::FRAC_PI_2 - atan(1.0 / x); }
    // For x > 0.4142 (tan(pi/8)), use atan(x) = pi/4 + atan((x-1)/(x+1))
    if x > 0.4142135623730950 {
        let reduced = (x - 1.0) / (x + 1.0);
        return constants::FRAC_PI_4 + atan_small(reduced);
    }
    atan_small(x)
}

/// Taylor series for atan, only for |x| < 0.42.
fn atan_small(x: f64) -> f64 {
    let x2 = x * x;
    let mut term = x;
    let mut sum = x;
    for n in 1..80 {
        term *= -x2;
        let contrib = term / (2 * n + 1) as f64;
        sum += contrib;
        if abs(contrib) < 1e-16 { break; }
    }
    sum
}

/// Two-argument arctangent.
pub fn atan2(y: f64, x: f64) -> f64 {
    if x > 0.0 { return atan(y / x); }
    if x < 0.0 && y >= 0.0 { return atan(y / x) + constants::PI; }
    if x < 0.0 && y < 0.0 { return atan(y / x) - constants::PI; }
    if x == 0.0 && y > 0.0 { return constants::FRAC_PI_2; }
    if x == 0.0 && y < 0.0 { return -constants::FRAC_PI_2; }
    0.0 // x == 0, y == 0
}

/// Inverse cotangent.
pub fn acot(x: f64) -> f64 { constants::FRAC_PI_2 - atan(x) }

/// Inverse secant.
pub fn asec(x: f64) -> f64 { acos(1.0 / x) }

/// Inverse cosecant.
pub fn acsc(x: f64) -> f64 { asin(1.0 / x) }

/// Degrees to radians.
#[inline]
pub fn to_radians(deg: f64) -> f64 { deg * constants::PI / 180.0 }

/// Radians to degrees.
#[inline]
pub fn to_degrees(rad: f64) -> f64 { rad * 180.0 / constants::PI }

/// Versine: 1 - cos(x).
pub fn versin(x: f64) -> f64 { 1.0 - cos(x) }

/// Haversine: (1 - cos(x)) / 2.
pub fn haversin(x: f64) -> f64 { (1.0 - cos(x)) / 2.0 }

/// Sinc: sin(x)/x, with sinc(0) = 1.
pub fn sinc(x: f64) -> f64 {
    if abs(x) < 1e-15 { 1.0 } else { sin(x) / x }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arithmetic::approx_eq;

    #[test]
    fn test_sin() {
        assert!(approx_eq(sin(0.0), 0.0, 1e-14));
        assert!(approx_eq(sin(constants::FRAC_PI_2), 1.0, 1e-14));
        assert!(approx_eq(sin(constants::PI), 0.0, 1e-14));
        assert!(approx_eq(sin(-constants::FRAC_PI_2), -1.0, 1e-14));
    }

    #[test]
    fn test_cos() {
        assert!(approx_eq(cos(0.0), 1.0, 1e-14));
        assert!(approx_eq(cos(constants::FRAC_PI_2), 0.0, 1e-14));
        assert!(approx_eq(cos(constants::PI), -1.0, 1e-14));
    }

    #[test]
    fn test_tan() {
        assert!(approx_eq(tan(0.0), 0.0, 1e-14));
        assert!(approx_eq(tan(constants::FRAC_PI_4), 1.0, 1e-12));
    }

    #[test]
    fn test_asin_acos() {
        assert!(approx_eq(asin(0.0), 0.0, 1e-14));
        assert!(approx_eq(asin(1.0), constants::FRAC_PI_2, 1e-12));
        assert!(approx_eq(acos(1.0), 0.0, 1e-14));
        assert!(approx_eq(acos(0.0), constants::FRAC_PI_2, 1e-12));
    }

    #[test]
    fn test_atan() {
        assert!(approx_eq(atan(0.0), 0.0, 1e-14));
        assert!(approx_eq(atan(1.0), constants::FRAC_PI_4, 1e-12));
    }

    #[test]
    fn test_atan2() {
        assert!(approx_eq(atan2(1.0, 1.0), constants::FRAC_PI_4, 1e-12));
        assert!(approx_eq(atan2(1.0, 0.0), constants::FRAC_PI_2, 1e-12));
    }

    #[test]
    fn test_sinc() {
        assert!(approx_eq(sinc(0.0), 1.0, 1e-14));
    }

    #[test]
    fn test_degrees_radians() {
        assert!(approx_eq(to_radians(180.0), constants::PI, 1e-14));
        assert!(approx_eq(to_degrees(constants::PI), 180.0, 1e-12));
    }
}
