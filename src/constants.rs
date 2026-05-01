//! Mathematical constants.

/// π — ratio of circle circumference to diameter
pub const PI: f64 = 3.14159265358979323846264338327950288;

/// τ = 2π — full circle in radians
pub const TAU: f64 = 6.28318530717958647692528676655900577;

/// e — Euler's number, base of natural logarithm
pub const E: f64 = 2.71828182845904523536028747135266250;

/// φ — Golden ratio (1 + √5) / 2
pub const PHI: f64 = 1.61803398874989484820458683436563812;

/// √2
pub const SQRT2: f64 = 1.41421356237309504880168872420969808;

/// √3
pub const SQRT3: f64 = 1.73205080756887729352744634150587237;

/// 1/√2
pub const FRAC_1_SQRT2: f64 = 0.70710678118654752440084436210484904;

/// π/2
pub const FRAC_PI_2: f64 = 1.57079632679489661923132169163975144;

/// π/3
pub const FRAC_PI_3: f64 = 1.04719755119659774615421446109316763;

/// π/4
pub const FRAC_PI_4: f64 = 0.78539816339744830961566084581987572;

/// π/6
pub const FRAC_PI_6: f64 = 0.52359877559829887307710723054658381;

/// 1/π
pub const FRAC_1_PI: f64 = 0.31830988618379067153776752674502872;

/// 2/π
pub const FRAC_2_PI: f64 = 0.63661977236758134307553505349005745;

/// 2/√π
pub const FRAC_2_SQRT_PI: f64 = 1.12837916709551257389615890312154517;

/// ln(2)
pub const LN2: f64 = 0.69314718055994530941723212145817657;

/// ln(10)
pub const LN10: f64 = 2.30258509299404568401799145468436421;

/// log₂(e)
pub const LOG2_E: f64 = 1.44269504088896340735992468100189214;

/// log₂(10)
pub const LOG2_10: f64 = 3.32192809488736234787031942948939018;

/// log₁₀(e)
pub const LOG10_E: f64 = 0.43429448190325182765112891891660508;

/// log₁₀(2)
pub const LOG10_2: f64 = 0.30102999566398119521373889472449302;

/// Euler–Mascheroni constant γ
pub const EULER_MASCHERONI: f64 = 0.57721566490153286060651209008240243;

/// Catalan's constant
pub const CATALAN: f64 = 0.91596559417721901505460351493238411;

/// Apéry's constant ζ(3)
pub const APERY: f64 = 1.20205690315959428539973816151144999;

/// Positive infinity
pub const INF: f64 = f64::INFINITY;

/// Negative infinity
pub const NEG_INF: f64 = f64::NEG_INFINITY;

/// Not a Number
pub const NAN: f64 = f64::NAN;

/// Machine epsilon for f64
pub const EPSILON: f64 = f64::EPSILON;

/// Maximum finite f64
pub const MAX: f64 = f64::MAX;

/// Minimum positive f64
pub const MIN_POSITIVE: f64 = f64::MIN_POSITIVE;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constants_relations() {
        assert!((TAU - 2.0 * PI).abs() < 1e-15);
        assert!((FRAC_PI_2 - PI / 2.0).abs() < 1e-15);
        assert!((FRAC_PI_4 - PI / 4.0).abs() < 1e-15);
        assert!((SQRT2 * FRAC_1_SQRT2 - 1.0).abs() < 1e-15);
        assert!((PHI * PHI - PHI - 1.0).abs() < 1e-14);
    }
}
