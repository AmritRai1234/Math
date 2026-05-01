//! Basic arithmetic and utility functions.

/// Absolute value.
#[inline]
pub fn abs(x: f64) -> f64 {
    if x < 0.0 { -x } else { x }
}

/// Sign of x: returns -1.0, 0.0, or 1.0.
#[inline]
pub fn sign(x: f64) -> f64 {
    if x > 0.0 {
        1.0
    } else if x < 0.0 {
        -1.0
    } else {
        0.0
    }
}

/// Copy the sign of `y` onto `x`.
#[inline]
pub fn copysign(x: f64, y: f64) -> f64 {
    abs(x) * sign(y)
}

/// Floor: largest integer ≤ x.
#[inline]
pub fn floor(x: f64) -> f64 {
    let i = x as i64;
    let f = i as f64;
    if x < f { f - 1.0 } else { f }
}

/// Ceiling: smallest integer ≥ x.
#[inline]
pub fn ceil(x: f64) -> f64 {
    let i = x as i64;
    let f = i as f64;
    if x > f { f + 1.0 } else { f }
}

/// Round to nearest integer (half rounds away from zero).
#[inline]
pub fn round(x: f64) -> f64 {
    if x >= 0.0 {
        floor(x + 0.5)
    } else {
        ceil(x - 0.5)
    }
}

/// Truncate toward zero.
#[inline]
pub fn trunc(x: f64) -> f64 {
    (x as i64) as f64
}

/// Fractional part: x - floor(x).
#[inline]
pub fn fract(x: f64) -> f64 {
    x - floor(x)
}

/// Clamp x to [min, max].
#[inline]
pub fn clamp(x: f64, min: f64, max: f64) -> f64 {
    if x < min {
        min
    } else if x > max {
        max
    } else {
        x
    }
}

/// Minimum of two values.
#[inline]
pub fn min(a: f64, b: f64) -> f64 {
    if a < b { a } else { b }
}

/// Maximum of two values.
#[inline]
pub fn max(a: f64, b: f64) -> f64 {
    if a > b { a } else { b }
}

/// Linear interpolation: a + t*(b - a).
#[inline]
pub fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + t * (b - a)
}

/// Inverse lerp: given v in [a,b], returns t in [0,1].
#[inline]
pub fn inverse_lerp(a: f64, b: f64, v: f64) -> f64 {
    if (b - a).abs() < f64::EPSILON {
        0.0
    } else {
        (v - a) / (b - a)
    }
}

/// Remap value from [in_min, in_max] to [out_min, out_max].
#[inline]
pub fn map_range(x: f64, in_min: f64, in_max: f64, out_min: f64, out_max: f64) -> f64 {
    let t = inverse_lerp(in_min, in_max, x);
    lerp(out_min, out_max, t)
}

/// Wrap x into [min, max) range.
#[inline]
pub fn wrap(x: f64, min_val: f64, max_val: f64) -> f64 {
    let range = max_val - min_val;
    if range.abs() < f64::EPSILON {
        return min_val;
    }
    let mut result = x - min_val;
    result = result - floor(result / range) * range;
    result + min_val
}

/// Step function: returns 0.0 if x < edge, else 1.0.
#[inline]
pub fn step(edge: f64, x: f64) -> f64 {
    if x < edge { 0.0 } else { 1.0 }
}

/// Smooth step (Hermite interpolation) between edge0 and edge1.
#[inline]
pub fn smoothstep(edge0: f64, edge1: f64, x: f64) -> f64 {
    let t = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// Ken Perlin's smoother step.
#[inline]
pub fn smootherstep(edge0: f64, edge1: f64, x: f64) -> f64 {
    let t = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
    t * t * t * (t * (t * 6.0 - 15.0) + 10.0)
}

/// Fused multiply-add: a * b + c.
#[inline]
pub fn fma(a: f64, b: f64, c: f64) -> f64 {
    a * b + c
}

/// Remainder (IEEE 754 style): x - round(x/y) * y.
#[inline]
pub fn remainder(x: f64, y: f64) -> f64 {
    x - round(x / y) * y
}

/// Modulo (always non-negative for positive divisor).
#[inline]
pub fn modulo(a: f64, b: f64) -> f64 {
    let r = a - floor(a / b) * b;
    if r < 0.0 { r + abs(b) } else { r }
}

/// Check if a value is NaN.
#[inline]
pub fn is_nan(x: f64) -> bool {
    x != x
}

/// Check if a value is infinite.
#[inline]
pub fn is_inf(x: f64) -> bool {
    x == f64::INFINITY || x == f64::NEG_INFINITY
}

/// Check if a value is finite.
#[inline]
pub fn is_finite(x: f64) -> bool {
    !is_nan(x) && !is_inf(x)
}

/// Check if two values are approximately equal within epsilon.
#[inline]
pub fn approx_eq(a: f64, b: f64, epsilon: f64) -> bool {
    abs(a - b) < epsilon
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_abs() {
        assert_eq!(abs(-5.0), 5.0);
        assert_eq!(abs(5.0), 5.0);
        assert_eq!(abs(0.0), 0.0);
    }

    #[test]
    fn test_sign() {
        assert_eq!(sign(42.0), 1.0);
        assert_eq!(sign(-42.0), -1.0);
        assert_eq!(sign(0.0), 0.0);
    }

    #[test]
    fn test_floor_ceil() {
        assert_eq!(floor(2.7), 2.0);
        assert_eq!(floor(-2.7), -3.0);
        assert_eq!(ceil(2.3), 3.0);
        assert_eq!(ceil(-2.3), -2.0);
    }

    #[test]
    fn test_round() {
        assert_eq!(round(2.5), 3.0);
        assert_eq!(round(2.4), 2.0);
        assert_eq!(round(-2.5), -3.0);
    }

    #[test]
    fn test_clamp() {
        assert_eq!(clamp(5.0, 0.0, 10.0), 5.0);
        assert_eq!(clamp(-1.0, 0.0, 10.0), 0.0);
        assert_eq!(clamp(15.0, 0.0, 10.0), 10.0);
    }

    #[test]
    fn test_lerp() {
        assert_eq!(lerp(0.0, 10.0, 0.5), 5.0);
        assert_eq!(lerp(0.0, 10.0, 0.0), 0.0);
        assert_eq!(lerp(0.0, 10.0, 1.0), 10.0);
    }

    #[test]
    fn test_map_range() {
        assert!(approx_eq(map_range(5.0, 0.0, 10.0, 0.0, 100.0), 50.0, 1e-10));
    }

    #[test]
    fn test_wrap() {
        assert!(approx_eq(wrap(370.0, 0.0, 360.0), 10.0, 1e-10));
        assert!(approx_eq(wrap(-10.0, 0.0, 360.0), 350.0, 1e-10));
    }

    #[test]
    fn test_smoothstep() {
        assert_eq!(smoothstep(0.0, 1.0, 0.0), 0.0);
        assert_eq!(smoothstep(0.0, 1.0, 1.0), 1.0);
        assert!(approx_eq(smoothstep(0.0, 1.0, 0.5), 0.5, 1e-10));
    }

    #[test]
    fn test_modulo() {
        assert!(approx_eq(modulo(7.0, 3.0), 1.0, 1e-10));
        assert!(approx_eq(modulo(-1.0, 3.0), 2.0, 1e-10));
    }
}
