//! Prediction models powered by differentiation and integration.
//!
//! This module provides forecasting tools that leverage numerical calculus:
//! - **Taylor extrapolation** — predict using derivative-based polynomial expansion
//! - **Derivative trend forecasting** — use rate-of-change and acceleration to project forward
//! - **Integral accumulation** — forecast cumulative quantities via numerical integration
//! - **ODE-based growth models** — exponential, logistic, and Gompertz growth via RK4
//! - **Lagrange extrapolation** — polynomial extrapolation from observed data
//! - **Ensemble prediction** — combine multiple methods with confidence intervals

use crate::arithmetic;
use crate::calculus;
use crate::interpolation;
use crate::regression;
use crate::powers;
use crate::exponential;

// ─── Derivative-Based Prediction ───────────────────────────────────────────

/// Estimate numerical derivatives up to order `max_order` from evenly-spaced data.
/// Returns a Vec where index i is the i-th derivative at the last data point.
pub fn estimate_derivatives(data: &[f64], dt: f64, max_order: usize) -> Vec<f64> {
    let mut derivs = Vec::with_capacity(max_order + 1);
    let mut current = data.to_vec();
    // 0th derivative = value itself (last point)
    derivs.push(*current.last().unwrap_or(&0.0));
    for _ in 0..max_order {
        if current.len() < 2 { derivs.push(0.0); continue; }
        let diff: Vec<f64> = current.windows(2).map(|w| (w[1] - w[0]) / dt).collect();
        derivs.push(*diff.last().unwrap_or(&0.0));
        current = diff;
    }
    derivs
}

/// Taylor series extrapolation: predict future values using derivatives at the last point.
///
/// Given data sampled at interval `dt`, estimates derivatives up to `order` and
/// uses the Taylor expansion: f(t₀+Δt) ≈ Σ f⁽ⁿ⁾(t₀) · Δtⁿ / n!
///
/// Returns `horizon` predicted values at intervals of `dt`.
pub fn taylor_predict(data: &[f64], dt: f64, order: usize, horizon: usize) -> Vec<f64> {
    let derivs = estimate_derivatives(data, dt, order);
    let mut result = Vec::with_capacity(horizon);
    for step in 1..=horizon {
        let delta = step as f64 * dt;
        let mut val = 0.0;
        let mut factorial = 1.0;
        for (n, &d) in derivs.iter().enumerate() {
            if n > 0 { factorial *= n as f64; }
            val += d * powers::pow(delta, n as f64) / factorial;
        }
        result.push(val);
    }
    result
}

/// Predict using instantaneous velocity and acceleration (2nd-order Taylor).
///
/// Simple physics-style prediction: x(t) = x₀ + v·t + ½·a·t²
pub fn kinematic_predict(data: &[f64], dt: f64, horizon: usize) -> Vec<f64> {
    taylor_predict(data, dt, 2, horizon)
}

// ─── Trend Slope Prediction ────────────────────────────────────────────────

/// Predict by fitting a linear trend to the last `window` points using their
/// numerical derivative (slope), then extrapolating forward.
pub fn slope_predict(data: &[f64], dt: f64, window: usize, horizon: usize) -> Vec<f64> {
    let w = window.min(data.len());
    let tail = &data[data.len() - w..];
    // Compute average slope over window
    let slopes: Vec<f64> = tail.windows(2).map(|p| (p[1] - p[0]) / dt).collect();
    let avg_slope = if slopes.is_empty() { 0.0 } else {
        slopes.iter().sum::<f64>() / slopes.len() as f64
    };
    let last = *data.last().unwrap_or(&0.0);
    (1..=horizon).map(|s| last + avg_slope * s as f64 * dt).collect()
}

/// Weighted slope prediction — recent slopes get exponentially higher weight.
pub fn weighted_slope_predict(data: &[f64], dt: f64, decay: f64, horizon: usize) -> Vec<f64> {
    if data.len() < 2 { return vec![*data.last().unwrap_or(&0.0); horizon]; }
    let slopes: Vec<f64> = data.windows(2).map(|w| (w[1] - w[0]) / dt).collect();
    let n = slopes.len();
    let mut wsum = 0.0;
    let mut wtotal = 0.0;
    for (i, &s) in slopes.iter().enumerate() {
        let w = exponential::exp(decay * (i as f64 - (n - 1) as f64));
        wsum += w * s;
        wtotal += w;
    }
    let slope = if wtotal > 0.0 { wsum / wtotal } else { 0.0 };
    let last = *data.last().unwrap();
    (1..=horizon).map(|s| last + slope * s as f64 * dt).collect()
}

// ─── Integration-Based Prediction ──────────────────────────────────────────

/// Predict cumulative total by integrating an estimated rate function.
///
/// Fits a polynomial of `degree` to the data (treated as a rate), then
/// integrates forward using Simpson's rule to get the cumulative forecast.
pub fn integral_predict(
    times: &[f64], rates: &[f64], degree: usize, t_start: f64, t_end: f64, steps: usize,
) -> f64 {
    let coeffs = regression::polynomial_regression(times, rates, degree);
    let rate_fn = |t: f64| -> f64 { regression::poly_predict(&coeffs, t) };
    calculus::integrate_simpson(&rate_fn, t_start, t_end, steps)
}

/// Forecast cumulative values at each future time step.
///
/// Given observed (time, rate) data, fits a polynomial rate model and integrates
/// to produce cumulative predictions at `horizon` future time steps of size `dt`.
pub fn cumulative_forecast(
    times: &[f64], rates: &[f64], degree: usize, dt: f64, horizon: usize,
) -> Vec<f64> {
    let coeffs = regression::polynomial_regression(times, rates, degree);
    let t0 = *times.last().unwrap_or(&0.0);
    let base: f64 = rates.iter().sum::<f64>(); // rough cumulative so far
    let rate_fn = |t: f64| -> f64 { regression::poly_predict(&coeffs, t) };
    let mut result = Vec::with_capacity(horizon);
    for step in 1..=horizon {
        let t_end = t0 + step as f64 * dt;
        let added = calculus::integrate_simpson(&rate_fn, t0, t_end, 100);
        result.push(base + added);
    }
    result
}

// ─── ODE Growth Models ─────────────────────────────────────────────────────

/// Exponential growth/decay prediction: dy/dt = r·y, solved with RK4.
///
/// - `y0`: initial value (last observed data point)
/// - `rate`: growth rate r (positive = growth, negative = decay)
/// - `dt`: time step
/// - `horizon`: number of future steps
pub fn exponential_growth(y0: f64, rate: f64, dt: f64, horizon: usize) -> Vec<f64> {
    let f = move |_t: f64, y: f64| -> f64 { rate * y };
    let result = crate::ode::rk4(&f, y0, 0.0, horizon as f64 * dt, dt);
    result.iter().skip(1).map(|&(_, y)| y).collect()
}

/// Logistic growth prediction: dy/dt = r·y·(1 - y/K), solved with RK4.
///
/// - `y0`: initial value
/// - `rate`: intrinsic growth rate
/// - `capacity`: carrying capacity K
pub fn logistic_growth(y0: f64, rate: f64, capacity: f64, dt: f64, horizon: usize) -> Vec<f64> {
    let f = move |_t: f64, y: f64| -> f64 { rate * y * (1.0 - y / capacity) };
    let result = crate::ode::rk4(&f, y0, 0.0, horizon as f64 * dt, dt);
    result.iter().skip(1).map(|&(_, y)| y).collect()
}

/// Gompertz growth: dy/dt = r·y·ln(K/y), solved with RK4.
pub fn gompertz_growth(y0: f64, rate: f64, capacity: f64, dt: f64, horizon: usize) -> Vec<f64> {
    let f = move |_t: f64, y: f64| -> f64 {
        if y <= 0.0 { return 0.0; }
        rate * y * exponential::ln(capacity / y)
    };
    let result = crate::ode::rk4(&f, y0, 0.0, horizon as f64 * dt, dt);
    result.iter().skip(1).map(|&(_, y)| y).collect()
}

/// Estimate growth rate from data using log-derivative: r ≈ d(ln y)/dt.
pub fn estimate_growth_rate(data: &[f64], dt: f64) -> f64 {
    let ln_data: Vec<f64> = data.iter()
        .filter(|&&v| v > 0.0)
        .map(|&v| exponential::ln(v))
        .collect();
    if ln_data.len() < 2 { return 0.0; }
    let slopes: Vec<f64> = ln_data.windows(2).map(|w| (w[1] - w[0]) / dt).collect();
    slopes.iter().sum::<f64>() / slopes.len() as f64
}

/// Auto-fit: detect whether data follows exponential or logistic growth,
/// estimate parameters, and return predictions.
pub fn auto_growth_predict(data: &[f64], dt: f64, horizon: usize) -> GrowthPrediction {
    let n = data.len();
    let y0 = *data.last().unwrap_or(&0.0);
    let rate = estimate_growth_rate(data, dt);

    // Check if growth is decelerating (logistic) by comparing first-half vs second-half rate
    let half = n / 2;
    let rate1 = if half >= 2 { estimate_growth_rate(&data[..half], dt) } else { rate };
    let rate2 = if n - half >= 2 { estimate_growth_rate(&data[half..], dt) } else { rate };

    if rate1 > 0.0 && rate2 > 0.0 && rate2 < rate1 * 0.7 {
        // Decelerating → logistic
        let capacity = y0 * 2.0 * (rate1 / (rate1 - rate2));
        let cap = if capacity > y0 { capacity } else { y0 * 10.0 };
        GrowthPrediction {
            model: GrowthModel::Logistic,
            rate,
            capacity: Some(cap),
            predictions: logistic_growth(y0, rate1, cap, dt, horizon),
        }
    } else {
        // Exponential
        GrowthPrediction {
            model: GrowthModel::Exponential,
            rate,
            capacity: None,
            predictions: exponential_growth(y0, rate, dt, horizon),
        }
    }
}

/// Growth model type detected by auto-fit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GrowthModel { Exponential, Logistic, Gompertz }

/// Result of auto growth prediction.
#[derive(Debug, Clone)]
pub struct GrowthPrediction {
    pub model: GrowthModel,
    pub rate: f64,
    pub capacity: Option<f64>,
    pub predictions: Vec<f64>,
}

// ─── Lagrange Extrapolation ────────────────────────────────────────────────

/// Extrapolate using Lagrange polynomial fitted to the last `window` data points.
pub fn lagrange_predict(data: &[f64], dt: f64, window: usize, horizon: usize) -> Vec<f64> {
    let w = window.min(data.len());
    let start = data.len() - w;
    let points: Vec<(f64, f64)> = (0..w)
        .map(|i| ((start + i) as f64 * dt, data[start + i]))
        .collect();
    let t_last = (data.len() - 1) as f64 * dt;
    (1..=horizon)
        .map(|s| interpolation::lagrange(&points, t_last + s as f64 * dt))
        .collect()
}

// ─── Ensemble Prediction ───────────────────────────────────────────────────

/// Combine multiple prediction methods and return mean prediction with
/// confidence interval (± 1 standard deviation across methods).
pub fn ensemble_predict(data: &[f64], dt: f64, horizon: usize) -> EnsembleForecast {
    let n = data.len();
    let window = n.min(8);
    // Collect predictions from different methods
    let p_taylor = taylor_predict(data, dt, 3, horizon);
    let p_slope = slope_predict(data, dt, window, horizon);
    let p_kinematic = kinematic_predict(data, dt, horizon);
    let p_weighted = weighted_slope_predict(data, dt, 2.0, horizon);

    let methods = [&p_taylor, &p_slope, &p_kinematic, &p_weighted];
    let mut mean = vec![0.0; horizon];
    let mut upper = vec![0.0; horizon];
    let mut lower = vec![0.0; horizon];

    for step in 0..horizon {
        let vals: Vec<f64> = methods.iter().map(|m| m[step]).collect();
        let m = vals.iter().sum::<f64>() / vals.len() as f64;
        let variance = vals.iter().map(|&v| (v - m) * (v - m)).sum::<f64>() / vals.len() as f64;
        let std = powers::sqrt(variance);
        mean[step] = m;
        upper[step] = m + std;
        lower[step] = m - std;
    }

    EnsembleForecast { mean, upper, lower }
}

/// Ensemble forecast result with confidence bands.
#[derive(Debug, Clone)]
pub struct EnsembleForecast {
    /// Mean prediction across all methods.
    pub mean: Vec<f64>,
    /// Upper confidence bound (mean + 1 std dev).
    pub upper: Vec<f64>,
    /// Lower confidence bound (mean - 1 std dev).
    pub lower: Vec<f64>,
}

// ─── Error Metrics ─────────────────────────────────────────────────────────

/// Mean Absolute Error between predictions and actual values.
pub fn mae(predicted: &[f64], actual: &[f64]) -> f64 {
    let n = predicted.len().min(actual.len());
    if n == 0 { return 0.0; }
    predicted.iter().zip(actual.iter())
        .map(|(&p, &a)| arithmetic::abs(p - a))
        .sum::<f64>() / n as f64
}

/// Root Mean Squared Error.
pub fn rmse(predicted: &[f64], actual: &[f64]) -> f64 {
    let n = predicted.len().min(actual.len());
    if n == 0 { return 0.0; }
    let mse: f64 = predicted.iter().zip(actual.iter())
        .map(|(&p, &a)| (p - a) * (p - a))
        .sum::<f64>() / n as f64;
    powers::sqrt(mse)
}

/// Mean Absolute Percentage Error.
pub fn mape(predicted: &[f64], actual: &[f64]) -> f64 {
    let n = predicted.len().min(actual.len());
    if n == 0 { return 0.0; }
    predicted.iter().zip(actual.iter())
        .filter(|(_, &a)| arithmetic::abs(a) > 1e-30)
        .map(|(&p, &a)| arithmetic::abs((a - p) / a))
        .sum::<f64>() / n as f64 * 100.0
}

// ─── Tests ─────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arithmetic::approx_eq;

    #[test]
    fn test_taylor_linear() {
        // Linear data: y = 2x → derivative = 2, should predict perfectly
        let data: Vec<f64> = (0..10).map(|i| 2.0 * i as f64).collect();
        let pred = taylor_predict(&data, 1.0, 1, 3);
        assert!(approx_eq(pred[0], 20.0, 0.1)); // next = 20
        assert!(approx_eq(pred[1], 22.0, 0.1)); // next = 22
        assert!(approx_eq(pred[2], 24.0, 0.1)); // next = 24
    }

    #[test]
    fn test_kinematic_quadratic() {
        // Quadratic: y = t² → v=2t, a=2 → predict y(5)=25, y(6)=36
        let data: Vec<f64> = (0..5).map(|i| (i as f64) * (i as f64)).collect();
        let pred = kinematic_predict(&data, 1.0, 2);
        // At t=4: y=16, v≈2*4=8, a≈2 → y(5)≈16+8+1=25
        assert!(approx_eq(pred[0], 25.0, 2.0));
    }

    #[test]
    fn test_slope_predict() {
        let data = [1.0, 3.0, 5.0, 7.0, 9.0];
        let pred = slope_predict(&data, 1.0, 5, 3);
        assert!(approx_eq(pred[0], 11.0, 0.01));
        assert!(approx_eq(pred[1], 13.0, 0.01));
    }

    #[test]
    fn test_exponential_growth() {
        let pred = exponential_growth(100.0, 0.1, 1.0, 5);
        // y(1) ≈ 100 * e^0.1 ≈ 110.517
        assert!(approx_eq(pred[0], 110.517, 0.1));
    }

    #[test]
    fn test_logistic_growth_bounded() {
        let pred = logistic_growth(10.0, 0.5, 100.0, 1.0, 50);
        // Should approach but not exceed capacity
        let last = *pred.last().unwrap();
        assert!(last < 100.0);
        assert!(last > 90.0); // should be near capacity after 50 steps
    }

    #[test]
    fn test_estimate_growth_rate() {
        // Pure exponential: y = e^(0.1*t)
        let data: Vec<f64> = (0..10).map(|i| exponential::exp(0.1 * i as f64)).collect();
        let r = estimate_growth_rate(&data, 1.0);
        assert!(approx_eq(r, 0.1, 0.01));
    }

    #[test]
    fn test_ensemble() {
        let data = [1.0, 2.0, 3.0, 4.0, 5.0];
        let forecast = ensemble_predict(&data, 1.0, 3);
        assert_eq!(forecast.mean.len(), 3);
        // Linear data → all methods should roughly agree → tight bands
        for i in 0..3 {
            assert!(forecast.upper[i] >= forecast.mean[i]);
            assert!(forecast.lower[i] <= forecast.mean[i]);
        }
    }

    #[test]
    fn test_mae_rmse() {
        let pred = [2.0, 4.0, 6.0];
        let actual = [2.5, 3.5, 6.5];
        assert!(approx_eq(mae(&pred, &actual), 0.5, 1e-10));
        assert!(approx_eq(rmse(&pred, &actual), 0.5, 0.01));
    }

    #[test]
    fn test_lagrange_predict() {
        // Quadratic: y = x²
        let data = [0.0, 1.0, 4.0, 9.0, 16.0];
        let pred = lagrange_predict(&data, 1.0, 5, 2);
        assert!(approx_eq(pred[0], 25.0, 0.1)); // 5² = 25
        assert!(approx_eq(pred[1], 36.0, 0.1)); // 6² = 36
    }

    #[test]
    fn test_integral_predict() {
        // Constant rate of 5.0 → integral from 0 to 10 = 50
        let times = [0.0, 1.0, 2.0, 3.0];
        let rates = [5.0, 5.0, 5.0, 5.0];
        let total = integral_predict(&times, &rates, 1, 0.0, 10.0, 100);
        assert!(approx_eq(total, 50.0, 0.1));
    }

    #[test]
    fn test_auto_growth_exponential() {
        let data: Vec<f64> = (0..10).map(|i| exponential::exp(0.1 * i as f64)).collect();
        let result = auto_growth_predict(&data, 1.0, 5);
        assert_eq!(result.model, GrowthModel::Exponential);
        assert!(result.rate > 0.0);
        assert_eq!(result.predictions.len(), 5);
    }

    #[test]
    fn test_weighted_slope() {
        let data = [1.0, 2.0, 3.0, 4.0, 5.0];
        let pred = weighted_slope_predict(&data, 1.0, 1.0, 3);
        // Constant slope of 1 → should predict 6, 7, 8
        assert!(approx_eq(pred[0], 6.0, 0.5));
    }
}
