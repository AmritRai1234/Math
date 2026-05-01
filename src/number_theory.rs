//! Number theory functions.

/// Greatest common divisor (Euclidean algorithm).
pub fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

/// Least common multiple.
pub fn lcm(a: u64, b: u64) -> u64 {
    if a == 0 || b == 0 { return 0; }
    a / gcd(a, b) * b
}

/// Extended Euclidean algorithm. Returns (gcd, x, y) such that a*x + b*y = gcd.
pub fn ext_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if a == 0 { return (b, 0, 1); }
    let (g, x, y) = ext_gcd(b % a, a);
    (g, y - (b / a) * x, x)
}

/// Check if two numbers are coprime.
pub fn is_coprime(a: u64, b: u64) -> bool {
    gcd(a, b) == 1
}

/// Modular exponentiation: (base^exp) mod modulus.
pub fn mod_pow(mut base: u64, mut exp: u64, modulus: u64) -> u64 {
    if modulus == 1 { return 0; }
    let mut result: u64 = 1;
    base %= modulus;
    while exp > 0 {
        if exp & 1 == 1 {
            result = mul_mod(result, base, modulus);
        }
        exp >>= 1;
        base = mul_mod(base, base, modulus);
    }
    result
}

/// Modular multiplication that avoids overflow.
fn mul_mod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

/// Modular inverse: a^(-1) mod m. Returns None if not invertible.
pub fn mod_inv(a: u64, m: u64) -> Option<u64> {
    let (g, x, _) = ext_gcd(a as i64, m as i64);
    if g != 1 { return None; }
    Some(((x % m as i64 + m as i64) % m as i64) as u64)
}

/// Primality test using Miller-Rabin.
pub fn is_prime(n: u64) -> bool {
    if n < 2 { return false; }
    if n == 2 || n == 3 || n == 5 || n == 7 { return true; }
    if n % 2 == 0 || n % 3 == 0 || n % 5 == 0 { return false; }

    // Write n-1 = d * 2^r
    let mut d = n - 1;
    let mut r = 0u32;
    while d % 2 == 0 { d /= 2; r += 1; }

    // Deterministic witnesses for n < 3,317,044,064,679,887,385,961,981
    let witnesses = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];

    'outer: for &a in &witnesses {
        if a >= n { continue; }
        let mut x = mod_pow(a, d, n);
        if x == 1 || x == n - 1 { continue; }
        for _ in 0..r - 1 {
            x = mul_mod(x, x, n);
            if x == n - 1 { continue 'outer; }
        }
        return false;
    }
    true
}

/// Sieve of Eratosthenes: returns all primes up to limit.
pub fn prime_sieve(limit: u64) -> Vec<u64> {
    if limit < 2 { return vec![]; }
    let n = limit as usize + 1;
    let mut sieve = vec![true; n];
    sieve[0] = false;
    if n > 1 { sieve[1] = false; }
    let mut i = 2;
    while i * i < n {
        if sieve[i] {
            let mut j = i * i;
            while j < n { sieve[j] = false; j += i; }
        }
        i += 1;
    }
    sieve.iter().enumerate()
        .filter(|&(_, &is_p)| is_p)
        .map(|(i, _)| i as u64)
        .collect()
}

/// Nth prime (1-indexed: nth_prime(1) = 2).
pub fn nth_prime(n: u64) -> u64 {
    if n == 0 { return 0; }
    let mut count = 0u64;
    let mut candidate = 1u64;
    while count < n {
        candidate += 1;
        if is_prime(candidate) { count += 1; }
    }
    candidate
}

/// Prime factorization: returns Vec of (prime, exponent).
pub fn prime_factors(mut n: u64) -> Vec<(u64, u32)> {
    let mut factors = Vec::new();
    if n <= 1 { return factors; }
    let mut d = 2u64;
    while d * d <= n {
        let mut count = 0u32;
        while n % d == 0 { n /= d; count += 1; }
        if count > 0 { factors.push((d, count)); }
        d += 1;
    }
    if n > 1 { factors.push((n, 1)); }
    factors
}

/// All divisors of n (unsorted).
pub fn divisors(n: u64) -> Vec<u64> {
    if n == 0 { return vec![]; }
    let mut divs = Vec::new();
    let mut i = 1u64;
    while i * i <= n {
        if n % i == 0 {
            divs.push(i);
            if i != n / i { divs.push(n / i); }
        }
        i += 1;
    }
    divs.sort();
    divs
}

/// Euler's totient function φ(n).
pub fn euler_totient(mut n: u64) -> u64 {
    if n == 0 { return 0; }
    let mut result = n;
    let mut d = 2u64;
    while d * d <= n {
        if n % d == 0 {
            while n % d == 0 { n /= d; }
            result -= result / d;
        }
        d += 1;
    }
    if n > 1 { result -= result / n; }
    result
}

/// Möbius function μ(n).
pub fn mobius(n: u64) -> i8 {
    if n == 1 { return 1; }
    let factors = prime_factors(n);
    for &(_, exp) in &factors {
        if exp > 1 { return 0; }
    }
    if factors.len() % 2 == 0 { 1 } else { -1 }
}

/// Check if n is a perfect number.
pub fn is_perfect(n: u64) -> bool {
    if n < 2 { return false; }
    let sum: u64 = divisors(n).iter().filter(|&&d| d != n).sum();
    sum == n
}

/// Sum of digits in base 10.
pub fn digit_sum(mut n: u64) -> u64 {
    let mut sum = 0;
    while n > 0 { sum += n % 10; n /= 10; }
    sum
}

/// Collatz sequence starting from n.
pub fn collatz(n: u64) -> Vec<u64> {
    let mut seq = vec![n];
    let mut current = n;
    while current != 1 && current != 0 {
        current = if current % 2 == 0 { current / 2 } else { 3 * current + 1 };
        seq.push(current);
    }
    seq
}

/// Jacobi symbol (a/n).
pub fn jacobi_symbol(mut a: i64, mut n: i64) -> i8 {
    if n <= 0 || n % 2 == 0 { return 0; }
    let mut result: i8 = 1;
    a = a.rem_euclid(n);
    while a != 0 {
        while a % 2 == 0 {
            a /= 2;
            let n_mod8 = n % 8;
            if n_mod8 == 3 || n_mod8 == 5 { result = -result; }
        }
        std::mem::swap(&mut a, &mut n);
        if a % 4 == 3 && n % 4 == 3 { result = -result; }
        a = a.rem_euclid(n);
    }
    if n == 1 { result } else { 0 }
}

/// Chinese Remainder Theorem: solve system x ≡ r[i] (mod m[i]).
pub fn chinese_remainder(remainders: &[u64], moduli: &[u64]) -> Option<u64> {
    if remainders.len() != moduli.len() || remainders.is_empty() { return None; }
    let mut total_mod: u64 = 1;
    for &m in moduli { total_mod = total_mod.checked_mul(m)?; }
    let mut result: u64 = 0;
    for i in 0..remainders.len() {
        let mi = total_mod / moduli[i];
        let inv = mod_inv(mi % moduli[i], moduli[i])?;
        result = (result + remainders[i] % total_mod * (mi % total_mod) % total_mod * (inv % total_mod)) % total_mod;
    }
    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gcd_lcm() {
        assert_eq!(gcd(12, 8), 4);
        assert_eq!(gcd(17, 13), 1);
        assert_eq!(lcm(4, 6), 12);
    }

    #[test]
    fn test_primes() {
        assert!(!is_prime(0));
        assert!(!is_prime(1));
        assert!(is_prime(2));
        assert!(is_prime(17));
        assert!(!is_prime(15));
        assert!(is_prime(104729));
    }

    #[test]
    fn test_sieve() {
        let primes = prime_sieve(20);
        assert_eq!(primes, vec![2, 3, 5, 7, 11, 13, 17, 19]);
    }

    #[test]
    fn test_mod_pow() {
        assert_eq!(mod_pow(2, 10, 1000), 24);
        assert_eq!(mod_pow(3, 4, 5), 1);
    }

    #[test]
    fn test_euler_totient() {
        assert_eq!(euler_totient(1), 1);
        assert_eq!(euler_totient(10), 4);
        assert_eq!(euler_totient(12), 4);
    }

    #[test]
    fn test_prime_factors() {
        assert_eq!(prime_factors(12), vec![(2, 2), (3, 1)]);
        assert_eq!(prime_factors(60), vec![(2, 2), (3, 1), (5, 1)]);
    }

    #[test]
    fn test_divisors() {
        assert_eq!(divisors(12), vec![1, 2, 3, 4, 6, 12]);
    }

    #[test]
    fn test_digit_sum() {
        assert_eq!(digit_sum(1234), 10);
    }

    #[test]
    fn test_is_perfect() {
        assert!(is_perfect(6));
        assert!(is_perfect(28));
        assert!(!is_perfect(12));
    }
}
