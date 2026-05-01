//! Power and root functions.

use crate::arithmetic::{abs, is_nan};
use crate::exponential;

/// Integer power: base^exp (exp is u32).
#[inline]
pub fn powi(mut base: f64, mut exp: u32) -> f64 {
    let mut result = 1.0;
    while exp > 0 {
        if exp & 1 == 1 { result *= base; }
        base *= base;
        exp >>= 1;
    }
    result
}

/// General power: base^exp for arbitrary f64 exponent.
pub fn pow(base: f64, exp: f64) -> f64 {
    if exp == 0.0 { return 1.0; }
    if base == 0.0 { return if exp > 0.0 { 0.0 } else { f64::INFINITY }; }
    if base == 1.0 { return 1.0; }
    let exp_int = exp as i64;
    if (exp - exp_int as f64).abs() < f64::EPSILON && exp_int >= 0 {
        return powi(base, exp_int as u32);
    }
    if base < 0.0 {
        if (exp - exp_int as f64).abs() > f64::EPSILON { return f64::NAN; }
        let result = exponential::exp(exp * exponential::ln(-base));
        if exp_int % 2 != 0 { -result } else { result }
    } else {
        exponential::exp(exp * exponential::ln(base))
    }
}

/// Square root using Newton-Raphson.
pub fn sqrt(x: f64) -> f64 {
    if x < 0.0 { return f64::NAN; }
    if x == 0.0 || is_nan(x) || x == f64::INFINITY { return x; }
    let i = x.to_bits();
    let i = 0x5fe6eb50c7b537a9 - (i >> 1);
    let mut guess = 1.0 / f64::from_bits(i);
    for _ in 0..8 { guess = 0.5 * (guess + x / guess); }
    guess
}

/// Fast inverse square root (1/√x).
pub fn fast_inv_sqrt(x: f64) -> f64 {
    if x <= 0.0 { return f64::NAN; }
    let mut guess = f64::from_bits(0x5fe6eb50c7b537a9 - (x.to_bits() >> 1));
    for _ in 0..4 { guess = guess * (1.5 - 0.5 * x * guess * guess); }
    guess
}

/// Cube root.
pub fn cbrt(x: f64) -> f64 {
    if x == 0.0 { return 0.0; }
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let ax = abs(x);
    let mut guess = if ax > 1.0 { ax / 3.0 } else { ax };
    for _ in 0..64 {
        let prev = guess;
        guess = (2.0 * guess + ax / (guess * guess)) / 3.0;
        if abs(guess - prev) < f64::EPSILON * abs(guess) { break; }
    }
    sign * guess
}

/// Nth root: x^(1/n).
pub fn nth_root(x: f64, n: u32) -> f64 {
    if n == 0 { return f64::NAN; }
    if n == 1 { return x; }
    if n == 2 { return sqrt(x); }
    if n == 3 { return cbrt(x); }
    if x == 0.0 { return 0.0; }
    if x < 0.0 && n % 2 == 0 { return f64::NAN; }
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let ax = abs(x);
    let nf = n as f64;
    let mut guess = exponential::exp(exponential::ln(ax) / nf);
    for _ in 0..64 {
        let prev = guess;
        let gn1 = powi(guess, n - 1);
        if gn1 == 0.0 { break; }
        guess = ((nf - 1.0) * guess + ax / gn1) / nf;
        if abs(guess - prev) < f64::EPSILON * abs(guess) { break; }
    }
    sign * guess
}

/// Hypotenuse: √(a² + b²) without overflow.
pub fn hypot(a: f64, b: f64) -> f64 {
    let a = abs(a);
    let b = abs(b);
    let (big, small) = if a > b { (a, b) } else { (b, a) };
    if big == 0.0 { return 0.0; }
    let ratio = small / big;
    big * sqrt(1.0 + ratio * ratio)
}

/// Integer square root: floor(√n).
pub fn isqrt(n: u64) -> u64 {
    if n == 0 { return 0; }
    let mut x = (n as f64).sqrt() as u64;
    while x > 0 && x * x > n { x -= 1; }
    while (x + 1) * (x + 1) <= n { x += 1; }
    x
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arithmetic::approx_eq;

    #[test]
    fn test_powi() {
        assert_eq!(powi(2.0, 10), 1024.0);
        assert_eq!(powi(3.0, 0), 1.0);
    }

    #[test]
    fn test_sqrt() {
        assert!(approx_eq(sqrt(4.0), 2.0, 1e-10));
        assert!(approx_eq(sqrt(2.0), 1.41421356237, 1e-10));
        assert!(sqrt(-1.0).is_nan());
    }

    #[test]
    fn test_cbrt() {
        assert!(approx_eq(cbrt(8.0), 2.0, 1e-10));
        assert!(approx_eq(cbrt(-8.0), -2.0, 1e-10));
    }

    #[test]
    fn test_hypot() {
        assert!(approx_eq(hypot(3.0, 4.0), 5.0, 1e-10));
    }

    #[test]
    fn test_isqrt() {
        assert_eq!(isqrt(9), 3);
        assert_eq!(isqrt(8), 2);
    }
}
