//! Probability distributions.

use crate::arithmetic::abs;
use crate::constants;
use crate::exponential;
use crate::powers;
use crate::special;

/// Uniform PDF on [a, b].
pub fn uniform_pdf(x: f64, a: f64, b: f64) -> f64 {
    if x >= a && x <= b { 1.0 / (b - a) } else { 0.0 }
}

/// Uniform CDF on [a, b].
pub fn uniform_cdf(x: f64, a: f64, b: f64) -> f64 {
    if x < a { 0.0 } else if x > b { 1.0 } else { (x - a) / (b - a) }
}

/// Standard normal (Gaussian) PDF: N(0, 1).
pub fn standard_normal_pdf(x: f64) -> f64 {
    let coeff = 1.0 / powers::sqrt(constants::TAU);
    coeff * exponential::exp(-0.5 * x * x)
}

/// Normal (Gaussian) PDF with mean μ and std dev σ.
pub fn normal_pdf(x: f64, mu: f64, sigma: f64) -> f64 {
    let z = (x - mu) / sigma;
    standard_normal_pdf(z) / sigma
}

/// Normal CDF using the error function.
pub fn normal_cdf(x: f64, mu: f64, sigma: f64) -> f64 {
    let z = (x - mu) / (sigma * constants::SQRT2);
    0.5 * (1.0 + special::erf(z))
}

/// Error function erf(x) — delegated to special module.
pub fn erf(x: f64) -> f64 {
    special::erf(x)
}

/// Complementary error function erfc(x).
pub fn erfc(x: f64) -> f64 {
    special::erfc(x)
}

/// Binomial PMF: P(X = k) for X ~ Bin(n, p).
pub fn binomial_pmf(k: u64, n: u64, p: f64) -> f64 {
    if k > n { return 0.0; }
    let coeff = crate::combinatorics::combinations(n, k) as f64;
    let pk = powers::pow(p, k as f64);
    let qnk = powers::pow(1.0 - p, (n - k) as f64);
    coeff * pk * qnk
}

/// Binomial CDF: P(X <= k) for X ~ Bin(n, p).
pub fn binomial_cdf(k: u64, n: u64, p: f64) -> f64 {
    let mut sum = 0.0;
    for i in 0..=k {
        sum += binomial_pmf(i, n, p);
    }
    sum
}

/// Poisson PMF: P(X = k) for X ~ Poisson(λ).
pub fn poisson_pmf(k: u64, lambda: f64) -> f64 {
    exponential::exp(-lambda) * powers::pow(lambda, k as f64) / crate::combinatorics::factorial(k) as f64
}

/// Exponential PDF: λ * e^(-λx) for x ≥ 0.
pub fn exponential_pdf(x: f64, lambda: f64) -> f64 {
    if x < 0.0 { 0.0 } else { lambda * exponential::exp(-lambda * x) }
}

/// Exponential CDF: 1 - e^(-λx) for x ≥ 0.
pub fn exponential_cdf(x: f64, lambda: f64) -> f64 {
    if x < 0.0 { 0.0 } else { 1.0 - exponential::exp(-lambda * x) }
}

/// Chi-squared PDF with k degrees of freedom.
pub fn chi_squared_pdf(x: f64, k: u64) -> f64 {
    if x < 0.0 { return 0.0; }
    let half_k = k as f64 / 2.0;
    let coeff = 1.0 / (powers::pow(2.0, half_k) * special::gamma(half_k));
    coeff * powers::pow(x, half_k - 1.0) * exponential::exp(-x / 2.0)
}

/// Student's t PDF with ν degrees of freedom.
pub fn student_t_pdf(x: f64, nu: f64) -> f64 {
    let coeff = special::gamma((nu + 1.0) / 2.0)
        / (powers::sqrt(nu * constants::PI) * special::gamma(nu / 2.0));
    coeff * powers::pow(1.0 + x * x / nu, -(nu + 1.0) / 2.0)
}

/// Beta PDF: x^(α-1) * (1-x)^(β-1) / B(α,β).
pub fn beta_pdf(x: f64, alpha: f64, beta: f64) -> f64 {
    if x < 0.0 || x > 1.0 { return 0.0; }
    let b = special::beta(alpha, beta);
    if abs(b) < 1e-300 { return f64::NAN; }
    powers::pow(x, alpha - 1.0) * powers::pow(1.0 - x, beta - 1.0) / b
}

/// Gamma PDF: x^(α-1) * e^(-x/β) / (β^α * Γ(α)).
pub fn gamma_pdf(x: f64, alpha: f64, beta_param: f64) -> f64 {
    if x < 0.0 { return 0.0; }
    let coeff = 1.0 / (powers::pow(beta_param, alpha) * special::gamma(alpha));
    coeff * powers::pow(x, alpha - 1.0) * exponential::exp(-x / beta_param)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arithmetic::approx_eq;

    #[test]
    fn test_standard_normal() {
        // PDF at 0 should be ~0.3989
        assert!(approx_eq(standard_normal_pdf(0.0), 0.39894228, 1e-6));
    }

    #[test]
    fn test_normal_cdf() {
        // CDF at mean should be 0.5
        assert!(approx_eq(normal_cdf(0.0, 0.0, 1.0), 0.5, 1e-6));
    }

    #[test]
    fn test_uniform() {
        assert!(approx_eq(uniform_pdf(0.5, 0.0, 1.0), 1.0, 1e-14));
        assert!(approx_eq(uniform_cdf(0.5, 0.0, 1.0), 0.5, 1e-14));
    }

    #[test]
    fn test_binomial() {
        // Coin flip: P(3 heads in 5 flips) = C(5,3) * 0.5^5 = 10/32
        assert!(approx_eq(binomial_pmf(3, 5, 0.5), 10.0 / 32.0, 1e-12));
    }

    #[test]
    fn test_poisson() {
        // λ=1, k=0: e^(-1) ≈ 0.3679
        assert!(approx_eq(poisson_pmf(0, 1.0), 0.36787944, 1e-6));
    }

    #[test]
    fn test_exponential() {
        assert!(approx_eq(exponential_pdf(0.0, 1.0), 1.0, 1e-14));
        assert!(approx_eq(exponential_cdf(0.0, 1.0), 0.0, 1e-14));
    }
}
