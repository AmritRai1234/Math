//! Hyperbolic functions.

use crate::arithmetic::abs;
use crate::exponential;

/// Hyperbolic sine: (e^x - e^(-x)) / 2.
pub fn sinh(x: f64) -> f64 {
    if abs(x) < 1e-5 { return x + x * x * x / 6.0; }
    let ex = exponential::exp(x);
    (ex - 1.0 / ex) / 2.0
}

/// Hyperbolic cosine: (e^x + e^(-x)) / 2.
pub fn cosh(x: f64) -> f64 {
    let ex = exponential::exp(x);
    (ex + 1.0 / ex) / 2.0
}

/// Hyperbolic tangent.
pub fn tanh(x: f64) -> f64 {
    if x > 20.0 { return 1.0; }
    if x < -20.0 { return -1.0; }
    let e2x = exponential::exp(2.0 * x);
    (e2x - 1.0) / (e2x + 1.0)
}

/// Hyperbolic cotangent.
pub fn coth(x: f64) -> f64 {
    if abs(x) < 1e-15 { return if x >= 0.0 { f64::INFINITY } else { f64::NEG_INFINITY }; }
    1.0 / tanh(x)
}

/// Hyperbolic secant.
pub fn sech(x: f64) -> f64 { 1.0 / cosh(x) }

/// Hyperbolic cosecant.
pub fn csch(x: f64) -> f64 { 1.0 / sinh(x) }

/// Inverse hyperbolic sine: ln(x + √(x²+1)).
pub fn asinh(x: f64) -> f64 {
    if abs(x) < 1e-8 { return x; }
    exponential::ln(x + crate::powers::sqrt(x * x + 1.0))
}

/// Inverse hyperbolic cosine: ln(x + √(x²-1)).
pub fn acosh(x: f64) -> f64 {
    if x < 1.0 { return f64::NAN; }
    exponential::ln(x + crate::powers::sqrt(x * x - 1.0))
}

/// Inverse hyperbolic tangent: 0.5 * ln((1+x)/(1-x)).
pub fn atanh(x: f64) -> f64 {
    if abs(x) >= 1.0 { return if x > 0.0 { f64::INFINITY } else { f64::NEG_INFINITY }; }
    0.5 * exponential::ln((1.0 + x) / (1.0 - x))
}

/// Inverse hyperbolic cotangent.
pub fn acoth(x: f64) -> f64 {
    if abs(x) <= 1.0 { return f64::NAN; }
    0.5 * exponential::ln((x + 1.0) / (x - 1.0))
}

/// Inverse hyperbolic secant.
pub fn asech(x: f64) -> f64 {
    if x <= 0.0 || x > 1.0 { return f64::NAN; }
    acosh(1.0 / x)
}

/// Inverse hyperbolic cosecant.
pub fn acsch(x: f64) -> f64 {
    if x == 0.0 { return f64::NAN; }
    asinh(1.0 / x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arithmetic::approx_eq;

    #[test]
    fn test_sinh_cosh() {
        assert!(approx_eq(sinh(0.0), 0.0, 1e-14));
        assert!(approx_eq(cosh(0.0), 1.0, 1e-14));
        // cosh² - sinh² = 1
        let x = 1.5;
        assert!(approx_eq(cosh(x) * cosh(x) - sinh(x) * sinh(x), 1.0, 1e-12));
    }

    #[test]
    fn test_tanh() {
        assert!(approx_eq(tanh(0.0), 0.0, 1e-14));
        assert!(tanh(100.0) < 1.001);
    }

    #[test]
    fn test_inverses() {
        let x = 2.0;
        assert!(approx_eq(sinh(asinh(x)), x, 1e-12));
        assert!(approx_eq(cosh(acosh(x)), x, 1e-12));
        let x = 0.5;
        assert!(approx_eq(tanh(atanh(x)), x, 1e-12));
    }
}
