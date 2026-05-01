//! Statistical functions.

use crate::arithmetic;
use crate::powers;

/// Sum of all elements.
pub fn sum(data: &[f64]) -> f64 {
    data.iter().fold(0.0, |acc, &x| acc + x)
}

/// Product of all elements.
pub fn product(data: &[f64]) -> f64 {
    data.iter().fold(1.0, |acc, &x| acc * x)
}

/// Arithmetic mean.
pub fn mean(data: &[f64]) -> f64 {
    if data.is_empty() { return f64::NAN; }
    sum(data) / data.len() as f64
}

/// Geometric mean.
pub fn geometric_mean(data: &[f64]) -> f64 {
    if data.is_empty() { return f64::NAN; }
    let n = data.len() as f64;
    let log_sum: f64 = data.iter().map(|&x| crate::exponential::ln(x)).sum();
    crate::exponential::exp(log_sum / n)
}

/// Harmonic mean.
pub fn harmonic_mean(data: &[f64]) -> f64 {
    if data.is_empty() { return f64::NAN; }
    let n = data.len() as f64;
    let inv_sum: f64 = data.iter().map(|&x| 1.0 / x).sum();
    n / inv_sum
}

/// Weighted mean.
pub fn weighted_mean(data: &[f64], weights: &[f64]) -> f64 {
    if data.len() != weights.len() || data.is_empty() { return f64::NAN; }
    let weighted: f64 = data.iter().zip(weights.iter()).map(|(&d, &w)| d * w).sum();
    let total_weight: f64 = weights.iter().sum();
    weighted / total_weight
}

/// Median (returns middle value of sorted data).
pub fn median(data: &[f64]) -> f64 {
    if data.is_empty() { return f64::NAN; }
    let mut sorted = data.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    if n % 2 == 0 {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    } else {
        sorted[n / 2]
    }
}

/// Mode (most frequent value). Returns first mode found.
pub fn mode(data: &[f64]) -> f64 {
    if data.is_empty() { return f64::NAN; }
    let mut sorted = data.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let mut best_val = sorted[0];
    let mut best_count = 1usize;
    let mut current_val = sorted[0];
    let mut current_count = 1usize;

    for &v in &sorted[1..] {
        if (v - current_val).abs() < f64::EPSILON {
            current_count += 1;
        } else {
            if current_count > best_count {
                best_count = current_count;
                best_val = current_val;
            }
            current_val = v;
            current_count = 1;
        }
    }
    if current_count > best_count { best_val = current_val; }
    best_val
}

/// Population variance.
pub fn variance(data: &[f64]) -> f64 {
    if data.is_empty() { return f64::NAN; }
    let m = mean(data);
    let sq_diff: f64 = data.iter().map(|&x| (x - m) * (x - m)).sum();
    sq_diff / data.len() as f64
}

/// Sample variance (Bessel's correction).
pub fn sample_variance(data: &[f64]) -> f64 {
    if data.len() < 2 { return f64::NAN; }
    let m = mean(data);
    let sq_diff: f64 = data.iter().map(|&x| (x - m) * (x - m)).sum();
    sq_diff / (data.len() - 1) as f64
}

/// Population standard deviation.
pub fn std_dev(data: &[f64]) -> f64 {
    powers::sqrt(variance(data))
}

/// Sample standard deviation.
pub fn sample_std_dev(data: &[f64]) -> f64 {
    powers::sqrt(sample_variance(data))
}

/// Covariance of two datasets.
pub fn covariance(x: &[f64], y: &[f64]) -> f64 {
    if x.len() != y.len() || x.is_empty() { return f64::NAN; }
    let mx = mean(x);
    let my = mean(y);
    let cov: f64 = x.iter().zip(y.iter()).map(|(&xi, &yi)| (xi - mx) * (yi - my)).sum();
    cov / x.len() as f64
}

/// Pearson correlation coefficient.
pub fn correlation(x: &[f64], y: &[f64]) -> f64 {
    let cov = covariance(x, y);
    let sx = std_dev(x);
    let sy = std_dev(y);
    if sx == 0.0 || sy == 0.0 { return f64::NAN; }
    cov / (sx * sy)
}

/// Skewness.
pub fn skewness(data: &[f64]) -> f64 {
    if data.len() < 3 { return f64::NAN; }
    let m = mean(data);
    let s = std_dev(data);
    let n = data.len() as f64;
    if s == 0.0 { return 0.0; }
    let m3: f64 = data.iter().map(|&x| ((x - m) / s).powi(3)).sum::<f64>();
    m3 / n
}

/// Kurtosis (excess kurtosis, 0 for normal distribution).
pub fn kurtosis(data: &[f64]) -> f64 {
    if data.len() < 4 { return f64::NAN; }
    let m = mean(data);
    let s = std_dev(data);
    let n = data.len() as f64;
    if s == 0.0 { return 0.0; }
    let m4: f64 = data.iter().map(|&x| ((x - m) / s).powi(4)).sum::<f64>();
    m4 / n - 3.0
}

/// Pth percentile (linear interpolation method).
pub fn percentile(data: &[f64], p: f64) -> f64 {
    if data.is_empty() || p < 0.0 || p > 100.0 { return f64::NAN; }
    let mut sorted = data.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let idx = (p / 100.0) * (sorted.len() - 1) as f64;
    let lo = idx as usize;
    let hi = (lo + 1).min(sorted.len() - 1);
    let frac = idx - lo as f64;
    sorted[lo] + frac * (sorted[hi] - sorted[lo])
}

/// Quartiles (Q1, Q2, Q3).
pub fn quartiles(data: &[f64]) -> (f64, f64, f64) {
    (percentile(data, 25.0), percentile(data, 50.0), percentile(data, 75.0))
}

/// Interquartile range.
pub fn iqr(data: &[f64]) -> f64 {
    let (q1, _, q3) = quartiles(data);
    q3 - q1
}

/// Z-score.
pub fn z_score(x: f64, mean: f64, std: f64) -> f64 {
    if std == 0.0 { return f64::NAN; }
    (x - mean) / std
}

/// Range: max - min.
pub fn range(data: &[f64]) -> f64 {
    if data.is_empty() { return f64::NAN; }
    let mn = data.iter().cloned().fold(f64::INFINITY, arithmetic::min);
    let mx = data.iter().cloned().fold(f64::NEG_INFINITY, arithmetic::max);
    mx - mn
}

/// Cumulative sum.
pub fn cumsum(data: &[f64]) -> Vec<f64> {
    let mut result = Vec::with_capacity(data.len());
    let mut acc = 0.0;
    for &x in data { acc += x; result.push(acc); }
    result
}

/// Moving average with window size w.
pub fn moving_avg(data: &[f64], w: usize) -> Vec<f64> {
    if w == 0 || data.len() < w { return vec![]; }
    let mut result = Vec::with_capacity(data.len() - w + 1);
    let mut window_sum: f64 = data[..w].iter().sum();
    result.push(window_sum / w as f64);
    for i in w..data.len() {
        window_sum += data[i] - data[i - w];
        result.push(window_sum / w as f64);
    }
    result
}

/// Minimum value in data.
pub fn min(data: &[f64]) -> f64 {
    if data.is_empty() { return f64::NAN; }
    data.iter().cloned().fold(f64::INFINITY, arithmetic::min)
}

/// Maximum value in data.
pub fn max(data: &[f64]) -> f64 {
    if data.is_empty() { return f64::NAN; }
    data.iter().cloned().fold(f64::NEG_INFINITY, arithmetic::max)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arithmetic::approx_eq;

    #[test]
    fn test_mean() {
        assert!(approx_eq(mean(&[1.0, 2.0, 3.0, 4.0, 5.0]), 3.0, 1e-14));
    }

    #[test]
    fn test_median() {
        assert!(approx_eq(median(&[1.0, 3.0, 2.0]), 2.0, 1e-14));
        assert!(approx_eq(median(&[1.0, 2.0, 3.0, 4.0]), 2.5, 1e-14));
    }

    #[test]
    fn test_variance() {
        let data = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        assert!(approx_eq(variance(&data), 4.0, 1e-12));
        assert!(approx_eq(std_dev(&data), 2.0, 1e-12));
    }

    #[test]
    fn test_correlation() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [2.0, 4.0, 6.0, 8.0, 10.0];
        assert!(approx_eq(correlation(&x, &y), 1.0, 1e-12));
    }

    #[test]
    fn test_cumsum() {
        assert_eq!(cumsum(&[1.0, 2.0, 3.0]), vec![1.0, 3.0, 6.0]);
    }

    #[test]
    fn test_moving_avg() {
        let r = moving_avg(&[1.0, 2.0, 3.0, 4.0, 5.0], 3);
        assert!(approx_eq(r[0], 2.0, 1e-14));
        assert!(approx_eq(r[1], 3.0, 1e-14));
        assert!(approx_eq(r[2], 4.0, 1e-14));
    }
}
