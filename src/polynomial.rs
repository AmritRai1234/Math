//! Polynomial arithmetic: evaluation, roots, differentiation, integration, multiplication.

use crate::arithmetic::abs;
use crate::powers;
use crate::complex::Complex;

/// A polynomial represented by coefficients [a0, a1, a2, ...] = a0 + a1*x + a2*x² + ...
#[derive(Clone, Debug, PartialEq)]
pub struct Polynomial {
    pub coeffs: Vec<f64>,
}

impl Polynomial {
    /// Create from coefficients [a0, a1, ..., an].
    pub fn new(coeffs: Vec<f64>) -> Self {
        let mut p = Self { coeffs };
        p.trim();
        p
    }

    /// Create from roots: (x - r0)(x - r1)...
    pub fn from_roots(roots: &[f64]) -> Self {
        let mut p = Self::new(vec![1.0]);
        for &r in roots {
            p = p.mul(&Self::new(vec![-r, 1.0]));
        }
        p
    }

    /// Zero polynomial.
    pub fn zero() -> Self { Self { coeffs: vec![0.0] } }

    /// Constant polynomial.
    pub fn constant(c: f64) -> Self { Self::new(vec![c]) }

    /// Degree of the polynomial.
    pub fn degree(&self) -> usize {
        if self.coeffs.is_empty() { 0 } else { self.coeffs.len().saturating_sub(1) }
    }

    /// Leading coefficient.
    pub fn leading_coeff(&self) -> f64 {
        *self.coeffs.last().unwrap_or(&0.0)
    }

    /// Remove trailing zero coefficients.
    fn trim(&mut self) {
        while self.coeffs.len() > 1 && abs(*self.coeffs.last().unwrap()) < 1e-30 {
            self.coeffs.pop();
        }
    }

    /// Evaluate polynomial at x using Horner's method.
    pub fn eval(&self, x: f64) -> f64 {
        let mut result = 0.0;
        for c in self.coeffs.iter().rev() {
            result = result * x + c;
        }
        result
    }

    /// Evaluate at complex x.
    pub fn eval_complex(&self, z: Complex) -> Complex {
        let mut result = Complex::new(0.0, 0.0);
        for c in self.coeffs.iter().rev() {
            result = result * z + Complex::new(*c, 0.0);
        }
        result
    }

    /// Formal derivative.
    pub fn derivative(&self) -> Self {
        if self.coeffs.len() <= 1 { return Self::zero(); }
        let coeffs: Vec<f64> = self.coeffs.iter().enumerate()
            .skip(1).map(|(i, &c)| c * i as f64).collect();
        Self::new(coeffs)
    }

    /// Nth derivative.
    pub fn nth_derivative(&self, n: u32) -> Self {
        let mut p = self.clone();
        for _ in 0..n { p = p.derivative(); }
        p
    }

    /// Formal integral with constant of integration c.
    pub fn integral(&self, c: f64) -> Self {
        let mut coeffs = vec![c];
        for (i, &coeff) in self.coeffs.iter().enumerate() {
            coeffs.push(coeff / (i + 1) as f64);
        }
        Self::new(coeffs)
    }

    /// Definite integral from a to b.
    pub fn definite_integral(&self, a: f64, b: f64) -> f64 {
        let anti = self.integral(0.0);
        anti.eval(b) - anti.eval(a)
    }

    /// Add two polynomials.
    pub fn add(&self, other: &Self) -> Self {
        let n = self.coeffs.len().max(other.coeffs.len());
        let mut coeffs = vec![0.0; n];
        for (i, &c) in self.coeffs.iter().enumerate() { coeffs[i] += c; }
        for (i, &c) in other.coeffs.iter().enumerate() { coeffs[i] += c; }
        Self::new(coeffs)
    }

    /// Subtract two polynomials.
    pub fn sub(&self, other: &Self) -> Self {
        let n = self.coeffs.len().max(other.coeffs.len());
        let mut coeffs = vec![0.0; n];
        for (i, &c) in self.coeffs.iter().enumerate() { coeffs[i] += c; }
        for (i, &c) in other.coeffs.iter().enumerate() { coeffs[i] -= c; }
        Self::new(coeffs)
    }

    /// Multiply two polynomials.
    pub fn mul(&self, other: &Self) -> Self {
        let n = self.coeffs.len() + other.coeffs.len() - 1;
        let mut coeffs = vec![0.0; n];
        for (i, &a) in self.coeffs.iter().enumerate() {
            for (j, &b) in other.coeffs.iter().enumerate() {
                coeffs[i + j] += a * b;
            }
        }
        Self::new(coeffs)
    }

    /// Scale by a constant.
    pub fn scale(&self, s: f64) -> Self {
        Self::new(self.coeffs.iter().map(|&c| c * s).collect())
    }

    /// Polynomial long division: returns (quotient, remainder).
    pub fn div_rem(&self, divisor: &Self) -> (Self, Self) {
        if divisor.leading_coeff() == 0.0 { panic!("Division by zero polynomial"); }
        let mut remainder = self.coeffs.clone();
        let dd = divisor.degree();
        let nd = self.degree();
        if nd < dd { return (Self::zero(), self.clone()); }
        let mut quotient = vec![0.0; nd - dd + 1];
        for i in (0..=nd - dd).rev() {
            quotient[i] = remainder[i + dd] / divisor.leading_coeff();
            for j in 0..=dd {
                remainder[i + j] -= quotient[i] * divisor.coeffs[j];
            }
        }
        (Self::new(quotient), Self::new(remainder))
    }

    /// Composition: self(other(x)).
    pub fn compose(&self, other: &Self) -> Self {
        let mut result = Self::zero();
        let mut power = Self::constant(1.0);
        for &c in &self.coeffs {
            result = result.add(&power.scale(c));
            power = power.mul(other);
        }
        result
    }

    /// Find real roots using Durand-Kerner method.
    pub fn roots(&self) -> Vec<Complex> {
        let n = self.degree();
        if n == 0 { return vec![]; }
        if n == 1 { return vec![Complex::new(-self.coeffs[0] / self.coeffs[1], 0.0)]; }
        if n == 2 {
            let a = self.coeffs[2]; let b = self.coeffs[1]; let c = self.coeffs[0];
            let disc = b * b - 4.0 * a * c;
            if disc >= 0.0 {
                let sd = powers::sqrt(disc);
                return vec![Complex::new((-b + sd) / (2.0 * a), 0.0), Complex::new((-b - sd) / (2.0 * a), 0.0)];
            } else {
                let sd = powers::sqrt(-disc);
                return vec![Complex::new(-b / (2.0 * a), sd / (2.0 * a)), Complex::new(-b / (2.0 * a), -sd / (2.0 * a))];
            }
        }

        // Durand-Kerner for higher degree
        let lc = self.leading_coeff();
        let monic: Vec<f64> = self.coeffs.iter().map(|&c| c / lc).collect();
        let monic_poly = Self::new(monic);

        let mut roots: Vec<Complex> = (0..n).map(|k| {
            let angle = 2.0 * crate::constants::PI * k as f64 / n as f64 + 0.4;
            Complex::new(0.4 * crate::trig::cos(angle), 0.4 * crate::trig::sin(angle))
        }).collect();

        for _ in 0..1000 {
            let mut max_delta = 0.0;
            for i in 0..n {
                let num = monic_poly.eval_complex(roots[i]);
                let mut denom = Complex::new(1.0, 0.0);
                for j in 0..n {
                    if i != j { denom = denom * (roots[i] - roots[j]); }
                }
                if denom.magnitude_sq() < 1e-60 { continue; }
                let delta = num / denom;
                roots[i] = roots[i] - delta;
                let d = delta.magnitude();
                if d > max_delta { max_delta = d; }
            }
            if max_delta < 1e-14 { break; }
        }
        roots
    }

    /// Find only real roots (filter complex roots with small imaginary part).
    pub fn real_roots(&self) -> Vec<f64> {
        self.roots().iter()
            .filter(|r| abs(r.im) < 1e-10)
            .map(|r| r.re)
            .collect()
    }
}

impl std::fmt::Display for Polynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut first = true;
        for (i, &c) in self.coeffs.iter().enumerate().rev() {
            if abs(c) < 1e-30 { continue; }
            if !first && c > 0.0 { write!(f, " + ")?; }
            else if !first && c < 0.0 { write!(f, " - ")?; }
            let ac = if first { c } else { abs(c) };
            match i {
                0 => write!(f, "{ac}")?,
                1 => if abs(ac - 1.0) < 1e-15 { write!(f, "x")? } else { write!(f, "{ac}x")? },
                _ => if abs(ac - 1.0) < 1e-15 { write!(f, "x^{i}")? } else { write!(f, "{ac}x^{i}")? },
            }
            first = false;
        }
        if first { write!(f, "0")?; }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arithmetic::approx_eq;

    #[test]
    fn test_eval() {
        let p = Polynomial::new(vec![1.0, 2.0, 3.0]); // 1 + 2x + 3x²
        assert!(approx_eq(p.eval(0.0), 1.0, 1e-14));
        assert!(approx_eq(p.eval(1.0), 6.0, 1e-14));
        assert!(approx_eq(p.eval(2.0), 17.0, 1e-14));
    }

    #[test]
    fn test_derivative() {
        let p = Polynomial::new(vec![1.0, 2.0, 3.0]); // 1 + 2x + 3x²
        let dp = p.derivative(); // 2 + 6x
        assert!(approx_eq(dp.eval(0.0), 2.0, 1e-14));
        assert!(approx_eq(dp.eval(1.0), 8.0, 1e-14));
    }

    #[test]
    fn test_integral() {
        let p = Polynomial::new(vec![2.0, 6.0]); // 2 + 6x
        assert!(approx_eq(p.definite_integral(0.0, 1.0), 5.0, 1e-14));
    }

    #[test]
    fn test_mul() {
        let a = Polynomial::new(vec![1.0, 1.0]); // 1 + x
        let b = Polynomial::new(vec![1.0, 1.0]); // 1 + x
        let c = a.mul(&b); // 1 + 2x + x²
        assert!(approx_eq(c.coeffs[0], 1.0, 1e-14));
        assert!(approx_eq(c.coeffs[1], 2.0, 1e-14));
        assert!(approx_eq(c.coeffs[2], 1.0, 1e-14));
    }

    #[test]
    fn test_roots_quadratic() {
        let p = Polynomial::new(vec![-6.0, 1.0, 1.0]); // x² + x - 6 = (x+3)(x-2)
        let mut roots = p.real_roots();
        roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!(approx_eq(roots[0], -3.0, 1e-10));
        assert!(approx_eq(roots[1], 2.0, 1e-10));
    }

    #[test]
    fn test_from_roots() {
        let p = Polynomial::from_roots(&[1.0, 2.0, 3.0]);
        assert!(approx_eq(p.eval(1.0), 0.0, 1e-10));
        assert!(approx_eq(p.eval(2.0), 0.0, 1e-10));
        assert!(approx_eq(p.eval(3.0), 0.0, 1e-10));
    }

    #[test]
    fn test_div_rem() {
        let a = Polynomial::new(vec![-4.0, 0.0, 1.0]); // x² - 4
        let b = Polynomial::new(vec![-2.0, 1.0]); // x - 2
        let (q, r) = a.div_rem(&b);
        assert!(approx_eq(q.eval(0.0), 2.0, 1e-10)); // quotient = x + 2
        assert!(approx_eq(r.eval(0.0), 0.0, 1e-10)); // remainder = 0
    }
}
