//! Exponential and logarithmic functions.

use crate::arithmetic::abs;
use crate::constants;

/// Exponential function e^x using Taylor series.
pub fn exp(x: f64) -> f64 {
    if x == 0.0 { return 1.0; }
    if x > 709.0 { return f64::INFINITY; }
    if x < -745.0 { return 0.0; }
    if x != x { return f64::NAN; }

    // Range reduction: e^x = 2^k * e^r where r = x - k*ln(2)
    let k = (x / constants::LN2).round() as i64;
    let r = x - (k as f64) * constants::LN2;

    // Taylor series for e^r (r is small)
    let mut sum = 1.0;
    let mut term = 1.0;
    for i in 1..30 {
        term *= r / (i as f64);
        sum += term;
        if abs(term) < f64::EPSILON * abs(sum) { break; }
    }

    // Multiply by 2^k using bit manipulation
    if k >= -1022 && k <= 1023 {
        let bits = ((k + 1023) as u64) << 52;
        sum * f64::from_bits(bits)
    } else {
        sum * pow2_f64(k)
    }
}

fn pow2_f64(k: i64) -> f64 {
    if k > 1023 { return f64::INFINITY; }
    if k < -1074 { return 0.0; }
    let mut result = 1.0;
    let mut remaining = k;
    while remaining > 0 { result *= 2.0; remaining -= 1; }
    while remaining < 0 { result *= 0.5; remaining += 1; }
    result
}

/// e^x - 1, accurate for small x.
pub fn expm1(x: f64) -> f64 {
    if abs(x) < 1e-5 {
        x + 0.5 * x * x + x * x * x / 6.0
    } else {
        exp(x) - 1.0
    }
}

/// 2^x
pub fn exp2(x: f64) -> f64 {
    exp(x * constants::LN2)
}

/// 10^x
pub fn exp10(x: f64) -> f64 {
    exp(x * constants::LN10)
}

/// Natural logarithm ln(x) using Halley's method.
pub fn ln(x: f64) -> f64 {
    if x < 0.0 { return f64::NAN; }
    if x == 0.0 { return f64::NEG_INFINITY; }
    if x == 1.0 { return 0.0; }
    if x == f64::INFINITY { return f64::INFINITY; }
    if x != x { return f64::NAN; }

    // Decompose x = m * 2^e where 1 <= m < 2
    let bits = x.to_bits();
    let e = ((bits >> 52) & 0x7FF) as i64 - 1023;
    let m_bits = (bits & 0x000FFFFFFFFFFFFF) | 0x3FF0000000000000;
    let m = f64::from_bits(m_bits);

    // ln(x) = e*ln(2) + ln(m)
    // For ln(m), use series: ln(m) = 2*sum((m-1)/(m+1))^(2k+1)/(2k+1))
    let s = (m - 1.0) / (m + 1.0);
    let s2 = s * s;
    let mut sum = s;
    let mut term = s;
    for k in 1..40 {
        term *= s2;
        let contrib = term / (2 * k + 1) as f64;
        sum += contrib;
        if abs(contrib) < f64::EPSILON * abs(sum) { break; }
    }
    (e as f64) * constants::LN2 + 2.0 * sum
}

/// ln(1 + x), accurate for small x.
pub fn ln1p(x: f64) -> f64 {
    if abs(x) < 1e-4 {
        x - 0.5 * x * x + x * x * x / 3.0 - x * x * x * x / 4.0
    } else {
        ln(1.0 + x)
    }
}

/// Log base 2.
pub fn log2(x: f64) -> f64 {
    ln(x) * constants::LOG2_E
}

/// Log base 10.
pub fn log10(x: f64) -> f64 {
    ln(x) * constants::LOG10_E
}

/// Log with arbitrary base.
pub fn log(x: f64, base: f64) -> f64 {
    ln(x) / ln(base)
}

/// Logistic sigmoid: 1 / (1 + e^(-x)).
pub fn logistic(x: f64) -> f64 {
    if x >= 0.0 {
        let e = exp(-x);
        1.0 / (1.0 + e)
    } else {
        let e = exp(x);
        e / (1.0 + e)
    }
}

/// Softplus: ln(1 + e^x).
pub fn softplus(x: f64) -> f64 {
    if x > 20.0 {
        x
    } else if x < -20.0 {
        0.0
    } else {
        ln(1.0 + exp(x))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arithmetic::approx_eq;

    #[test]
    fn test_exp() {
        assert!(approx_eq(exp(0.0), 1.0, 1e-14));
        assert!(approx_eq(exp(1.0), constants::E, 1e-14));
        assert!(approx_eq(exp(2.0), 7.38905609893065, 1e-10));
        assert!(approx_eq(exp(-1.0), 1.0 / constants::E, 1e-14));
    }

    #[test]
    fn test_ln() {
        assert!(approx_eq(ln(1.0), 0.0, 1e-14));
        assert!(approx_eq(ln(constants::E), 1.0, 1e-14));
        assert!(approx_eq(ln(2.0), constants::LN2, 1e-14));
        assert!(approx_eq(ln(10.0), constants::LN10, 1e-12));
    }

    #[test]
    fn test_log2() {
        assert!(approx_eq(log2(1.0), 0.0, 1e-14));
        assert!(approx_eq(log2(2.0), 1.0, 1e-12));
        assert!(approx_eq(log2(8.0), 3.0, 1e-12));
    }

    #[test]
    fn test_log10() {
        assert!(approx_eq(log10(1.0), 0.0, 1e-14));
        assert!(approx_eq(log10(10.0), 1.0, 1e-12));
        assert!(approx_eq(log10(100.0), 2.0, 1e-12));
    }

    #[test]
    fn test_logistic() {
        assert!(approx_eq(logistic(0.0), 0.5, 1e-14));
        assert!(logistic(100.0) < 1.001);
        assert!(logistic(-100.0) > -0.001);
    }
}
