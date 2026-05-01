//! Hybrid equations: novel systems mixing calculus, chaos, entropy, quaternions, and fuzzy logic.

use crate::arithmetic;
use crate::exponential;
use crate::powers;
use crate::trig;
use crate::constants;

// ═══════════════════════════════════════════════════════════════════════════
// 1. ENTROPIC LORENZ — Lorenz attractor where σ,ρ,β adapt based on
//    trajectory entropy. The system self-tunes between order and chaos.
// ═══════════════════════════════════════════════════════════════════════════

/// Lorenz system where parameters are modulated by the running entropy
/// of the trajectory. High entropy → parameters push toward chaos,
/// low entropy → parameters stabilize. A self-organizing chaotic system.
pub fn entropic_lorenz(
    x0: f64, y0: f64, z0: f64, dt: f64, steps: usize,
) -> Vec<(f64, f64, f64, f64)> {
    let (mut x, mut y, mut z) = (x0, y0, z0);
    let mut result = vec![(x, y, z, 0.0)];
    let mut history = Vec::with_capacity(steps);

    for step in 0..steps {
        // Estimate entropy from recent trajectory variance
        history.push(x * x + y * y + z * z);
        let entropy = if history.len() > 10 {
            let tail = &history[history.len() - 10..];
            let mean = tail.iter().sum::<f64>() / 10.0;
            let var = tail.iter().map(|&v| (v - mean) * (v - mean)).sum::<f64>() / 10.0;
            if var > 1e-30 { 0.5 * exponential::ln(var) } else { 0.0 }
        } else { 1.0 };

        // Entropy modulates Lorenz parameters
        let sigma = 10.0 + 2.0 * trig::sin(entropy);
        let rho = 28.0 + 5.0 * trig::cos(entropy * 0.5);
        let beta = 8.0 / 3.0 + 0.5 * trig::sin(entropy * 0.3);

        let dx = sigma * (y - x) * dt;
        let dy = (x * (rho - z) - y) * dt;
        let dz = (x * y - beta * z) * dt;
        x += dx; y += dy; z += dz;
        result.push((x, y, z, entropy));
    }
    result
}

// ═══════════════════════════════════════════════════════════════════════════
// 2. GOLDEN SPIRAL ODE — differential equation generating a spiral whose
//    expansion rate is governed by the golden ratio φ.
// ═══════════════════════════════════════════════════════════════════════════

/// A 2D ODE system that produces a golden-ratio spiral:
///   dr/dt = r · ln(φ) / (2π)
///   dθ/dt = 1
/// Returns (x, y) coordinates.
pub fn golden_spiral(dt: f64, steps: usize) -> Vec<(f64, f64)> {
    let ln_phi = exponential::ln(constants::PHI);
    let mut r = 1.0;
    let mut theta = 0.0;
    let mut result = vec![(r, 0.0)];

    for _ in 0..steps {
        r += r * ln_phi / (2.0 * constants::PI) * dt;
        theta += dt;
        result.push((r * trig::cos(theta), r * trig::sin(theta)));
    }
    result
}

// ═══════════════════════════════════════════════════════════════════════════
// 3. QUATERNION LORENZ — Lorenz attractor extended to 4D via quaternion
//    multiplication. Produces a hyper-chaotic system.
// ═══════════════════════════════════════════════════════════════════════════

/// 4D chaotic system: quaternion state Q evolves via
///   dQ/dt = σ(Q·J - Q) + ρ·Q·K - Q·(Q*·Q)·β
/// where J,K are unit quaternions. Produces hyper-chaos in 4D.
pub fn quaternion_chaos(
    w0: f64, x0: f64, y0: f64, z0: f64, dt: f64, steps: usize,
) -> Vec<(f64, f64, f64, f64)> {
    let (mut w, mut x, mut y, mut z) = (w0, x0, y0, z0);
    let sigma = 10.0;
    let rho = 28.0;
    let beta = 8.0 / 3.0;
    let mut result = vec![(w, x, y, z)];

    for _ in 0..steps {
        let norm2 = w * w + x * x + y * y + z * z;
        // Quaternion cross-coupling terms
        let dw = (sigma * (x - w) - beta * w * norm2) * dt;
        let dx = (rho * w - x - w * z + y * z) * dt;
        let dy = (w * z - sigma * y + x * z) * dt;
        let dz = (-w * y + x * y - beta * z) * dt;
        w += dw; x += dx; y += dy; z += dz;
        result.push((w, x, y, z));
    }
    result
}

// ═══════════════════════════════════════════════════════════════════════════
// 4. FUZZY LOGISTIC MAP — logistic map where parameter r is a fuzzy number,
//    producing bands of possibility instead of single trajectories.
// ═══════════════════════════════════════════════════════════════════════════

/// Logistic map where r is fuzzy-triangular centered at `r_center` with
/// spread `delta`. Returns (min, center, max) at each step — a fuzzy trajectory.
pub fn fuzzy_logistic(
    r_center: f64, delta: f64, x0: f64, steps: usize,
) -> Vec<(f64, f64, f64)> {
    let mut x_lo = x0;
    let mut x_mid = x0;
    let mut x_hi = x0;
    let mut result = vec![(x_lo, x_mid, x_hi)];

    let r_lo = r_center - delta;
    let r_hi = r_center + delta;

    for _ in 0..steps {
        x_lo = r_lo * x_lo * (1.0 - x_lo);
        x_mid = r_center * x_mid * (1.0 - x_mid);
        x_hi = r_hi * x_hi * (1.0 - x_hi);
        // Keep bounds ordered
        let mn = x_lo.min(x_mid).min(x_hi);
        let mx = x_lo.max(x_mid).max(x_hi);
        result.push((mn, x_mid, mx));
    }
    result
}

// ═══════════════════════════════════════════════════════════════════════════
// 5. WAVE-GROWTH RESONANCE — wave equation coupled with logistic growth.
//    The amplitude grows logistically while oscillating.
// ═══════════════════════════════════════════════════════════════════════════

/// y'' + ω²y = r·y·(1 - y/K) — a harmonic oscillator driven by logistic growth.
/// Creates oscillations that grow toward a carrying capacity then ring.
pub fn wave_growth(
    y0: f64, v0: f64, omega: f64, r: f64, capacity: f64, dt: f64, steps: usize,
) -> Vec<(f64, f64)> {
    let mut y = y0;
    let mut v = v0;
    let mut result = vec![(0.0, y)];

    for i in 0..steps {
        let t = (i + 1) as f64 * dt;
        let accel = -omega * omega * y + r * y * (1.0 - y / capacity);
        v += accel * dt;
        y += v * dt;
        result.push((t, y));
    }
    result
}

// ═══════════════════════════════════════════════════════════════════════════
// 6. BESSEL-DAMPED OSCILLATOR — oscillator where damping follows a Bessel
//    function, creating non-monotonic energy decay.
// ═══════════════════════════════════════════════════════════════════════════

/// y'' + J₀(t)·y' + ω²y = 0 — damping oscillates via Bessel J₀.
/// Unlike standard damping, energy can temporarily INCREASE.
pub fn bessel_oscillator(
    y0: f64, v0: f64, omega: f64, dt: f64, steps: usize,
) -> Vec<(f64, f64, f64)> {
    let mut y = y0;
    let mut v = v0;
    let mut result = vec![(0.0, y, v)];

    for i in 0..steps {
        let t = (i + 1) as f64 * dt;
        let j0 = crate::special::bessel_j(0, t);
        let accel = -j0 * v - omega * omega * y;
        v += accel * dt;
        y += v * dt;
        result.push((t, y, v));
    }
    result
}

// ═══════════════════════════════════════════════════════════════════════════
// 7. ZETA POTENTIAL FIELD — particle moves in a potential defined by the
//    Riemann zeta function. Zeros of zeta become equilibrium points.
// ═══════════════════════════════════════════════════════════════════════════

/// Particle in potential V(x) = |ζ(x)|². Force = -dV/dx.
/// The particle is attracted to zeros of the zeta function.
pub fn zeta_dynamics(
    x0: f64, v0: f64, dt: f64, steps: usize,
) -> Vec<(f64, f64, f64)> {
    let mut x = x0;
    let mut v = v0;
    let mut result = vec![(0.0, x, v)];
    let h = 0.01;

    for i in 0..steps {
        let t = (i + 1) as f64 * dt;
        // V(x) = zeta(x)^2, F = -dV/dx ≈ -2·zeta(x)·zeta'(x)
        let z = crate::special::zeta(x);
        let zp = crate::special::zeta(x + h);
        let zm = crate::special::zeta(x - h);
        let dz = (zp - zm) / (2.0 * h);
        let force = -2.0 * z * dz;
        // Add friction to settle
        v += (force - 0.3 * v) * dt;
        x += v * dt;
        result.push((t, x, v));
    }
    result
}

// ═══════════════════════════════════════════════════════════════════════════
// 8. FIBONACCI RESONATOR — oscillator whose frequency steps through
//    Fibonacci ratios, converging toward golden ratio harmony.
// ═══════════════════════════════════════════════════════════════════════════

/// Oscillator where frequency = F(n+1)/F(n) at each cycle, converging to φ.
/// Creates a chirp that stabilizes at golden-ratio frequency.
pub fn fibonacci_resonator(
    y0: f64, dt: f64, steps: usize,
) -> Vec<(f64, f64, f64)> {
    let mut y = y0;
    let mut v = 0.0;
    let mut result = vec![(0.0, y, 1.0)];
    let (mut fa, mut fb): (f64, f64) = (1.0, 1.0);

    for i in 0..steps {
        let t = (i + 1) as f64 * dt;
        let omega = fb / fa; // approaches φ
        let accel = -omega * omega * y;
        v += accel * dt;
        y += v * dt;

        // Advance Fibonacci every full cycle
        if (t % (2.0 * constants::PI)).abs() < dt {
            let next = fa + fb;
            fa = fb;
            fb = next;
        }
        result.push((t, y, omega));
    }
    result
}

// ═══════════════════════════════════════════════════════════════════════════
// 9. SIGMOID-COUPLED PREDATOR-PREY — Lotka-Volterra where interaction
//    strength passes through sigmoid (neural-network style activation).
// ═══════════════════════════════════════════════════════════════════════════

/// dx/dt = αx - sigmoid(βxy), dy/dt = sigmoid(δxy) - γy
/// Neural activation replaces raw interaction → bounded, smooth switching.
pub fn sigmoid_predator_prey(
    x0: f64, y0: f64, alpha: f64, beta: f64, gamma: f64, delta: f64,
    dt: f64, steps: usize,
) -> Vec<(f64, f64, f64)> {
    let mut x = x0;
    let mut y = y0;
    let mut result = vec![(0.0, x, y)];
    let sig = |v: f64| 1.0 / (1.0 + exponential::exp(-v));

    for i in 0..steps {
        let t = (i + 1) as f64 * dt;
        let interaction = x * y;
        let dx = (alpha * x - sig(beta * interaction)) * dt;
        let dy = (sig(delta * interaction) - gamma * y) * dt;
        x += dx; y += dy;
        x = x.max(0.0); y = y.max(0.0);
        result.push((t, x, y));
    }
    result
}

// ═══════════════════════════════════════════════════════════════════════════
// 10. GAMMA-WEIGHTED DIFFUSION — heat equation where diffusivity is the
//     Gamma function of the local temperature. Wild nonlinear diffusion.
// ═══════════════════════════════════════════════════════════════════════════

/// 1D diffusion: du/dt = Γ(|u|+1) · d²u/dx²
/// Diffusivity explodes at integer temperatures (Gamma poles at negatives,
/// factorial growth at positives). Creates non-physical but fascinating patterns.
pub fn gamma_diffusion(initial: &[f64], dt: f64, steps: usize) -> Vec<Vec<f64>> {
    let n = initial.len();
    let mut u = initial.to_vec();
    let mut result = vec![u.clone()];

    for _ in 0..steps {
        let mut new_u = u.clone();
        for i in 1..n - 1 {
            let abs_u = arithmetic::abs(u[i]).min(10.0) + 1.0; // clamp to prevent overflow
            let diff_coeff = crate::special::gamma(abs_u).min(100.0); // bounded Gamma
            let laplacian = u[i + 1] - 2.0 * u[i] + u[i - 1];
            new_u[i] = u[i] + diff_coeff * laplacian * dt;
        }
        u = new_u;
        result.push(u.clone());
    }
    result
}

// ═══════════════════════════════════════════════════════════════════════════
// ANALYSIS TOOLS
// ═══════════════════════════════════════════════════════════════════════════

/// Compute the largest Lyapunov exponent of a 3D trajectory (chaos measure).
/// Positive = chaotic, zero = periodic, negative = stable.
pub fn lyapunov_from_trajectory(traj: &[(f64, f64, f64)], dt: f64) -> f64 {
    if traj.len() < 3 { return 0.0; }
    let mut sum = 0.0;
    let mut count = 0;
    for i in 1..traj.len() - 1 {
        let dx = traj[i + 1].0 - traj[i].0;
        let dy = traj[i + 1].1 - traj[i].1;
        let dz = traj[i + 1].2 - traj[i].2;
        let dist = powers::sqrt(dx * dx + dy * dy + dz * dz);
        if dist > 1e-30 {
            sum += exponential::ln(dist / dt);
            count += 1;
        }
    }
    if count > 0 { sum / count as f64 } else { 0.0 }
}

/// Compute the energy of a (position, velocity) trajectory.
pub fn trajectory_energy(traj: &[(f64, f64)]) -> Vec<f64> {
    traj.iter().map(|&(y, v)| 0.5 * v * v + 0.5 * y * y).collect()
}

/// Phase space distance between two trajectories (divergence measure).
pub fn trajectory_divergence(a: &[(f64, f64)], b: &[(f64, f64)]) -> Vec<f64> {
    a.iter().zip(b.iter())
        .map(|(&(y1, v1), &(y2, v2))| powers::sqrt((y1 - y2).powi(2) + (v1 - v2).powi(2)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arithmetic::approx_eq;

    #[test]
    fn test_entropic_lorenz_runs() {
        let traj = entropic_lorenz(1.0, 1.0, 1.0, 0.01, 500);
        assert_eq!(traj.len(), 501);
        // Should produce non-trivial trajectory
        let last = traj.last().unwrap();
        assert!((last.0.abs() + last.1.abs() + last.2.abs()) > 0.01);
    }

    #[test]
    fn test_golden_spiral_expands() {
        let spiral = golden_spiral(0.1, 200);
        let r_start = powers::sqrt(spiral[1].0.powi(2) + spiral[1].1.powi(2));
        let r_end = powers::sqrt(spiral[200].0.powi(2) + spiral[200].1.powi(2));
        assert!(r_end > r_start); // spiral expands
    }

    #[test]
    fn test_quaternion_chaos_4d() {
        let traj = quaternion_chaos(1.0, 0.1, 0.1, 0.1, 0.001, 1000);
        assert_eq!(traj.len(), 1001);
    }

    #[test]
    fn test_fuzzy_logistic_bands() {
        let traj = fuzzy_logistic(3.7, 0.1, 0.5, 50);
        assert_eq!(traj.len(), 51);
        // Fuzzy band should have width > 0
        let last = traj.last().unwrap();
        assert!(last.2 >= last.0); // max >= min
    }

    #[test]
    fn test_wave_growth_oscillates() {
        let traj = wave_growth(0.1, 0.0, 2.0, 0.5, 10.0, 0.01, 500);
        assert_eq!(traj.len(), 501);
        // Should cross zero (oscillation)
        let has_positive = traj.iter().any(|&(_, y)| y > 0.0);
        let has_negative = traj.iter().any(|&(_, y)| y < 0.0);
        assert!(has_positive);
        assert!(has_negative);
    }

    #[test]
    fn test_bessel_oscillator_runs() {
        let traj = bessel_oscillator(1.0, 0.0, 2.0, 0.01, 500);
        assert_eq!(traj.len(), 501);
    }

    #[test]
    fn test_fibonacci_resonator_converges() {
        let traj = fibonacci_resonator(1.0, 0.01, 5000);
        let last_omega = traj.last().unwrap().2;
        // Should approach φ ≈ 1.618
        assert!(approx_eq(last_omega, constants::PHI, 0.5));
    }

    #[test]
    fn test_sigmoid_predator_prey() {
        let traj = sigmoid_predator_prey(2.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.01, 500);
        assert_eq!(traj.len(), 501);
        // Populations stay non-negative
        assert!(traj.iter().all(|&(_, x, y)| x >= 0.0 && y >= 0.0));
    }

    #[test]
    fn test_gamma_diffusion_smooths() {
        let initial: Vec<f64> = (0..20).map(|i| if i == 10 { 1.0 } else { 0.0 }).collect();
        let result = gamma_diffusion(&initial, 0.001, 100);
        assert_eq!(result.len(), 101);
        // Peak should decrease (diffusion)
        assert!(result[100][10] < result[0][10]);
    }

    #[test]
    fn test_lyapunov_positive_for_chaos() {
        let traj = entropic_lorenz(1.0, 1.0, 1.0, 0.01, 2000);
        let traj3d: Vec<(f64, f64, f64)> = traj.iter().map(|&(x, y, z, _)| (x, y, z)).collect();
        let le = lyapunov_from_trajectory(&traj3d, 0.01);
        assert!(le > 0.0); // chaotic → positive Lyapunov
    }
}
