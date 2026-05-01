//! Special mathematical functions: gamma, beta, Bessel, zeta, error function, etc.

use crate::arithmetic::abs;
use crate::constants;
use crate::exponential;
use crate::powers;
use crate::trig;

// ─── Gamma & Related ───

/// Gamma function Γ(x) using Lanczos approximation.
pub fn gamma(x: f64) -> f64 {
    if x != x { return f64::NAN; }
    if x <= 0.0 {
        let frac = x - (x as i64 as f64);
        if abs(frac) < 1e-15 { return f64::INFINITY; } // poles at non-positive integers
        // Reflection: Γ(x) = π / (sin(πx) * Γ(1-x))
        return constants::PI / (trig::sin(constants::PI * x) * gamma(1.0 - x));
    }

    // Lanczos approximation with g=7
    const P: [f64; 9] = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ];

    let x = x - 1.0;
    let mut sum = P[0];
    for i in 1..9 {
        sum += P[i] / (x + i as f64);
    }

    let t = x + 7.5;
    constants::SQRT2 * powers::sqrt(constants::PI) * powers::pow(t, x + 0.5) * exponential::exp(-t) * sum
}

/// Log-gamma: ln(|Γ(x)|).
pub fn lgamma(x: f64) -> f64 {
    if x <= 0.0 { return f64::INFINITY; }
    exponential::ln(abs(gamma(x)))
}

/// Digamma function ψ(x) = d/dx ln(Γ(x)).
pub fn digamma(x: f64) -> f64 {
    if x <= 0.0 && abs(x - (x as i64 as f64)) < 1e-15 { return f64::NAN; }

    let mut result = 0.0;
    let mut z = x;

    // Recurrence to shift to large z
    while z < 8.0 {
        result -= 1.0 / z;
        z += 1.0;
    }

    // Asymptotic expansion
    let z2 = 1.0 / (z * z);
    result += exponential::ln(z) - 0.5 / z
        - z2 * (1.0 / 12.0 - z2 * (1.0 / 120.0 - z2 * (1.0 / 252.0)));
    result
}

/// Beta function B(a, b) = Γ(a)Γ(b) / Γ(a+b).
pub fn beta(a: f64, b: f64) -> f64 {
    exponential::exp(lgamma(a) + lgamma(b) - lgamma(a + b))
}

// ─── Error Function ───

/// Error function erf(x) using Horner form approximation.
pub fn erf(x: f64) -> f64 {
    if x == 0.0 { return 0.0; }
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = abs(x);

    // Abramowitz & Stegun approximation 7.1.26
    let p = 0.3275911;
    let a1 = 0.254829592;
    let a2 = -0.284496736;
    let a3 = 1.421413741;
    let a4 = -1.453152027;
    let a5 = 1.061405429;

    let t = 1.0 / (1.0 + p * x);
    let poly = t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * a5))));
    sign * (1.0 - poly * exponential::exp(-x * x))
}

/// Complementary error function erfc(x) = 1 - erf(x).
pub fn erfc(x: f64) -> f64 {
    1.0 - erf(x)
}

/// Inverse error function (approximation).
pub fn erf_inv(x: f64) -> f64 {
    if x <= -1.0 { return f64::NEG_INFINITY; }
    if x >= 1.0 { return f64::INFINITY; }

    let a = 0.147;
    let ln_term = exponential::ln(1.0 - x * x);
    let part1 = 2.0 / (constants::PI * a) + ln_term / 2.0;
    let part2 = ln_term / a;

    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    sign * powers::sqrt(powers::sqrt(part1 * part1 - part2) - part1)
}

// ─── Incomplete Gamma & Beta ───

/// Lower incomplete gamma γ(a, x) using series expansion.
pub fn inc_gamma_lower(a: f64, x: f64) -> f64 {
    if x < 0.0 { return f64::NAN; }
    if x == 0.0 { return 0.0; }

    let mut sum = 0.0;
    let mut term = 1.0 / a;
    sum += term;
    for n in 1..200 {
        term *= x / (a + n as f64);
        sum += term;
        if abs(term) < 1e-15 * abs(sum) { break; }
    }
    sum * powers::pow(x, a) * exponential::exp(-x)
}

/// Regularized incomplete gamma P(a, x) = γ(a,x) / Γ(a).
pub fn inc_gamma(a: f64, x: f64) -> f64 {
    inc_gamma_lower(a, x) / gamma(a)
}

/// Regularized incomplete beta I_x(a, b) — simple series approximation.
pub fn inc_beta(x: f64, a: f64, b: f64) -> f64 {
    if x < 0.0 || x > 1.0 { return f64::NAN; }
    if x == 0.0 { return 0.0; }
    if x == 1.0 { return 1.0; }

    // Use symmetry for x > 0.5
    if x > (a + 1.0) / (a + b + 2.0) {
        return 1.0 - inc_beta(1.0 - x, b, a);
    }

    let lbeta = lgamma(a) + lgamma(b) - lgamma(a + b);
    let front = exponential::exp(a * exponential::ln(x) + b * exponential::ln(1.0 - x) - lbeta) / a;

    // Continued fraction (Lentz's method simplified)
    let mut f = 1.0;
    let mut c = 1.0;
    let mut d;
    for m in 1..200 {
        let mf = m as f64;
        // Even step
        let num_even = mf * (b - mf) * x / ((a + 2.0 * mf - 1.0) * (a + 2.0 * mf));
        d = 1.0 + num_even / 1.0;
        c = 1.0 + num_even / c;
        if abs(d) < 1e-30 { d = 1e-30; }
        if abs(c) < 1e-30 { c = 1e-30; }
        f *= c / d;

        // Odd step
        let num_odd = -((a + mf) * (a + b + mf) * x) / ((a + 2.0 * mf) * (a + 2.0 * mf + 1.0));
        d = 1.0 + num_odd / 1.0;
        c = 1.0 + num_odd / c;
        if abs(d) < 1e-30 { d = 1e-30; }
        if abs(c) < 1e-30 { c = 1e-30; }
        let delta = c / d;
        f *= delta;
        if abs(delta - 1.0) < 1e-12 { break; }
    }

    front * f
}

// ─── Bessel Functions ───

/// Bessel function of the first kind J_n(x) using series expansion.
pub fn bessel_j(n: i32, x: f64) -> f64 {
    let nf = n.unsigned_abs() as f64;
    let mut sum = 0.0;
    let mut term = powers::pow(x / 2.0, nf) / gamma(nf + 1.0);
    sum += term;
    for k in 1..60 {
        term *= -(x * x) / (4.0 * k as f64 * (nf + k as f64));
        sum += term;
        if abs(term) < 1e-16 * abs(sum) { break; }
    }
    if n < 0 && n % 2 != 0 { -sum } else { sum }
}

/// Bessel function of the second kind Y_n(x).
pub fn bessel_y(n: i32, x: f64) -> f64 {
    if x <= 0.0 { return f64::NAN; }
    // Y_n = (J_n(x)*cos(nπ) - J_{-n}(x)) / sin(nπ)
    // For integer n, use limit form
    let nf = n as f64;
    let cos_npi = trig::cos(nf * constants::PI);
    let sin_npi = trig::sin(nf * constants::PI);
    if abs(sin_npi) > 1e-10 {
        (bessel_j(n, x) * cos_npi - bessel_j(-n, x)) / sin_npi
    } else {
        // For integer n, use numerical differentiation approximation
        let eps = 1e-6;
        let j_plus = bessel_j(n, x + eps);
        let j_minus = bessel_j(n, x - eps);
        let dj = (j_plus - j_minus) / (2.0 * eps);
        // Rough Neumann function approximation
        (2.0 / constants::PI) * (dj * exponential::ln(x / 2.0) + bessel_j(n, x) / x)
    }
}

/// Modified Bessel function of the first kind I_n(x).
pub fn bessel_i(n: i32, x: f64) -> f64 {
    let nf = n.unsigned_abs() as f64;
    let mut sum = 0.0;
    let mut term = powers::pow(x / 2.0, nf) / gamma(nf + 1.0);
    sum += term;
    for k in 1..60 {
        term *= (x * x) / (4.0 * k as f64 * (nf + k as f64));
        sum += term;
        if abs(term) < 1e-16 * abs(sum) { break; }
    }
    sum
}

/// Modified Bessel function of the second kind K_n(x) (approximate).
pub fn bessel_k(n: i32, x: f64) -> f64 {
    if x <= 0.0 { return f64::INFINITY; }
    let nf = n as f64;
    constants::FRAC_PI_2 * (bessel_i(-n.abs(), x) - bessel_i(n.abs(), x))
        / trig::sin(nf.abs() * constants::PI + 0.0001) // avoid exact zero
}

// ─── Riemann Zeta ───

/// Riemann zeta function ζ(s) for real s > 1 (and approximate elsewhere).
pub fn zeta(s: f64) -> f64 {
    if s == 1.0 { return f64::INFINITY; }

    if s > 1.0 {
        // Direct summation with Euler-Maclaurin
        let n = 100;
        let mut sum = 0.0;
        for k in 1..=n {
            sum += 1.0 / powers::pow(k as f64, s);
        }
        // Integral tail correction: ∫_n^∞ x^(-s) dx = n^(1-s)/(s-1)
        sum += powers::pow(n as f64, 1.0 - s) / (s - 1.0);
        sum += 0.5 * powers::pow(n as f64, -s); // Euler-Maclaurin correction
        return sum;
    }

    // For s < 1, use functional equation: ζ(s) = 2^s π^(s-1) sin(πs/2) Γ(1-s) ζ(1-s)
    if s < 0.0 {
        let one_minus_s = 1.0 - s;
        return powers::pow(2.0, s)
            * powers::pow(constants::PI, s - 1.0)
            * trig::sin(constants::PI * s / 2.0)
            * gamma(one_minus_s)
            * zeta(one_minus_s);
    }

    // 0 < s < 1: Dirichlet eta function η(s) = (1-2^(1-s))ζ(s)
    let mut eta = 0.0;
    for k in 1..200 {
        let sign = if k % 2 == 1 { 1.0 } else { -1.0 };
        eta += sign / powers::pow(k as f64, s);
    }
    eta / (1.0 - powers::pow(2.0, 1.0 - s))
}

/// Hurwitz zeta function ζ(s, a) = Σ 1/(n+a)^s.
pub fn hurwitz_zeta(s: f64, a: f64) -> f64 {
    if s <= 1.0 || a <= 0.0 { return f64::NAN; }
    let mut sum = 0.0;
    for n in 0..1000 {
        let term = 1.0 / powers::pow(n as f64 + a, s);
        sum += term;
        if abs(term) < 1e-15 * abs(sum) { break; }
    }
    sum
}

// ─── Polylogarithm ───

/// Polylogarithm Li_s(z) = Σ z^k / k^s.
pub fn polylog(s: f64, z: f64) -> f64 {
    if abs(z) > 1.0 { return f64::NAN; } // series only converges for |z| <= 1
    let mut sum = 0.0;
    for k in 1..500 {
        let term = powers::pow(z, k as f64) / powers::pow(k as f64, s);
        sum += term;
        if abs(term) < 1e-15 * abs(sum) { break; }
    }
    sum
}

// ─── Lambert W ───

/// Lambert W function (principal branch W₀): W(x) * e^W(x) = x.
pub fn lambert_w(x: f64) -> f64 {
    if x < -1.0 / constants::E { return f64::NAN; }
    if abs(x) < 1e-10 { return x; }

    // Initial guess
    let mut w = if x < 1.0 {
        x
    } else {
        exponential::ln(x) - exponential::ln(exponential::ln(x) + 1.0)
    };

    // Halley's method
    for _ in 0..50 {
        let ew = exponential::exp(w);
        let wew = w * ew;
        let f = wew - x;
        let fp = ew * (w + 1.0);
        let fpp = ew * (w + 2.0);
        let dw = f / (fp - f * fpp / (2.0 * fp));
        w -= dw;
        if abs(dw) < 1e-15 * abs(w) { break; }
    }
    w
}

// ─── Airy Functions ───

/// Airy Ai function (series approximation for small x).
pub fn airy_ai(x: f64) -> f64 {
    // Ai(x) = c1 * f(x) - c2 * g(x)
    // where f and g are power series
    let c1 = 1.0 / (powers::pow(3.0, 2.0 / 3.0) * gamma(2.0 / 3.0));
    let c2 = 1.0 / (powers::pow(3.0, 1.0 / 3.0) * gamma(1.0 / 3.0));

    let mut f = 1.0;
    let mut g = x;
    let mut f_term = 1.0;
    let mut g_term = x;
    let x3 = x * x * x;

    for k in 1..30 {
        f_term *= x3 / ((3 * k) as f64 * (3 * k - 1) as f64);
        g_term *= x3 / ((3 * k + 1) as f64 * (3 * k) as f64);
        f += f_term;
        g += g_term;
        if abs(f_term) + abs(g_term) < 1e-16 { break; }
    }

    c1 * f - c2 * g
}

/// Airy Bi function (series approximation for small x).
pub fn airy_bi(x: f64) -> f64 {
    let c1 = 1.0 / (powers::pow(3.0, 1.0 / 6.0) * gamma(2.0 / 3.0));
    let c2 = powers::pow(3.0, 1.0 / 6.0) / gamma(1.0 / 3.0);

    let mut f = 1.0;
    let mut g = x;
    let mut f_term = 1.0;
    let mut g_term = x;
    let x3 = x * x * x;

    for k in 1..30 {
        f_term *= x3 / ((3 * k) as f64 * (3 * k - 1) as f64);
        g_term *= x3 / ((3 * k + 1) as f64 * (3 * k) as f64);
        f += f_term;
        g += g_term;
        if abs(f_term) + abs(g_term) < 1e-16 { break; }
    }

    c1 * f + c2 * g
}

// ─── Orthogonal Polynomials ───

/// Legendre polynomial P_n(x) using recurrence.
pub fn legendre_p(n: u32, x: f64) -> f64 {
    if n == 0 { return 1.0; }
    if n == 1 { return x; }
    let mut p0 = 1.0;
    let mut p1 = x;
    for k in 2..=n {
        let kf = k as f64;
        let p2 = ((2.0 * kf - 1.0) * x * p1 - (kf - 1.0) * p0) / kf;
        p0 = p1;
        p1 = p2;
    }
    p1
}

/// Chebyshev polynomial T_n(x) using recurrence.
pub fn chebyshev_t(n: u32, x: f64) -> f64 {
    if n == 0 { return 1.0; }
    if n == 1 { return x; }
    let mut t0 = 1.0;
    let mut t1 = x;
    for _ in 2..=n {
        let t2 = 2.0 * x * t1 - t0;
        t0 = t1;
        t1 = t2;
    }
    t1
}

/// Hermite polynomial H_n(x) (physicist's convention).
pub fn hermite(n: u32, x: f64) -> f64 {
    if n == 0 { return 1.0; }
    if n == 1 { return 2.0 * x; }
    let mut h0 = 1.0;
    let mut h1 = 2.0 * x;
    for k in 2..=n {
        let h2 = 2.0 * x * h1 - 2.0 * (k - 1) as f64 * h0;
        h0 = h1;
        h1 = h2;
    }
    h1
}

/// Laguerre polynomial L_n(x) using recurrence.
pub fn laguerre(n: u32, x: f64) -> f64 {
    if n == 0 { return 1.0; }
    if n == 1 { return 1.0 - x; }
    let mut l0 = 1.0;
    let mut l1 = 1.0 - x;
    for k in 2..=n {
        let kf = k as f64;
        let l2 = ((2.0 * kf - 1.0 - x) * l1 - (kf - 1.0) * l0) / kf;
        l0 = l1;
        l1 = l2;
    }
    l1
}

// ─── Elliptic Integrals ───

/// Complete elliptic integral of the first kind K(m).
pub fn elliptic_k(m: f64) -> f64 {
    if m >= 1.0 { return f64::INFINITY; }
    if m < 0.0 { return f64::NAN; }
    // AGM method
    let mut a = 1.0;
    let mut b = powers::sqrt(1.0 - m);
    for _ in 0..50 {
        let a_new = (a + b) / 2.0;
        let b_new = powers::sqrt(a * b);
        if abs(a_new - b_new) < 1e-15 { a = a_new; break; }
        a = a_new;
        b = b_new;
    }
    constants::FRAC_PI_2 / a
}

/// Complete elliptic integral of the second kind E(m).
pub fn elliptic_e(m: f64) -> f64 {
    if m < 0.0 { return f64::NAN; }
    if m == 0.0 { return constants::FRAC_PI_2; }
    if m == 1.0 { return 1.0; }
    // Series expansion
    let mut sum = 1.0;
    let mut term = 1.0;
    for n in 1..50 {
        let nf = n as f64;
        term *= (2.0 * nf - 1.0) / (2.0 * nf) * (2.0 * nf - 1.0) / (2.0 * nf);
        sum -= term * powers::pow(m, nf) / (2.0 * nf - 1.0);
        if abs(term) < 1e-16 { break; }
    }
    sum * constants::FRAC_PI_2
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arithmetic::approx_eq;

    #[test]
    fn test_gamma() {
        assert!(approx_eq(gamma(1.0), 1.0, 1e-12));
        assert!(approx_eq(gamma(5.0), 24.0, 1e-10)); // 4!
        assert!(approx_eq(gamma(0.5), powers::sqrt(constants::PI), 1e-10));
    }

    #[test]
    fn test_beta() {
        // B(1,1) = 1
        assert!(approx_eq(beta(1.0, 1.0), 1.0, 1e-10));
    }

    #[test]
    fn test_erf() {
        assert!(approx_eq(erf(0.0), 0.0, 1e-14));
        assert!(erf(3.0) > 0.999);
        assert!(erf(-3.0) < -0.999);
    }

    #[test]
    fn test_zeta() {
        assert!(approx_eq(zeta(2.0), constants::PI * constants::PI / 6.0, 1e-4));
    }

    #[test]
    fn test_legendre() {
        assert!(approx_eq(legendre_p(0, 0.5), 1.0, 1e-14));
        assert!(approx_eq(legendre_p(1, 0.5), 0.5, 1e-14));
        assert!(approx_eq(legendre_p(2, 0.5), -0.125, 1e-14));
    }

    #[test]
    fn test_chebyshev() {
        assert!(approx_eq(chebyshev_t(0, 0.5), 1.0, 1e-14));
        assert!(approx_eq(chebyshev_t(1, 0.5), 0.5, 1e-14));
    }
}
