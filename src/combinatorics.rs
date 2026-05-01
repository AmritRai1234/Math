//! Combinatorics functions.

/// Factorial n!
pub fn factorial(n: u64) -> u64 {
    if n <= 1 { return 1; }
    let mut result: u64 = 1;
    for i in 2..=n {
        result = result.saturating_mul(i);
    }
    result
}

/// Double factorial n!!
pub fn double_factorial(n: u64) -> u64 {
    if n <= 1 { return 1; }
    let mut result: u64 = 1;
    let mut i = n;
    while i > 1 {
        result = result.saturating_mul(i);
        i -= 2;
    }
    result
}

/// Falling factorial: n * (n-1) * ... * (n-k+1).
pub fn falling_factorial(n: u64, k: u64) -> u64 {
    let mut result: u64 = 1;
    for i in 0..k {
        result = result.saturating_mul(n - i);
    }
    result
}

/// Rising factorial (Pochhammer): n * (n+1) * ... * (n+k-1).
pub fn rising_factorial(n: u64, k: u64) -> u64 {
    let mut result: u64 = 1;
    for i in 0..k {
        result = result.saturating_mul(n + i);
    }
    result
}

/// Permutations P(n, k) = n! / (n-k)!
pub fn permutations(n: u64, k: u64) -> u64 {
    if k > n { return 0; }
    falling_factorial(n, k)
}

/// Combinations C(n, k) = n! / (k! * (n-k)!)
pub fn combinations(n: u64, k: u64) -> u64 {
    if k > n { return 0; }
    let k = if k > n - k { n - k } else { k };
    let mut result: u64 = 1;
    for i in 0..k {
        result = result / (i + 1) * (n - i) + result % (i + 1) * (n - i) / (i + 1);
    }
    result
}

/// Binomial coefficient (alias for combinations).
pub fn binomial_coeff(n: u64, k: u64) -> u64 {
    combinations(n, k)
}

/// Multinomial coefficient: n! / (k1! * k2! * ... * km!)
pub fn multinomial(n: u64, ks: &[u64]) -> u64 {
    let sum: u64 = ks.iter().sum();
    if sum != n { return 0; }
    let mut result = factorial(n);
    for &k in ks {
        result /= factorial(k);
    }
    result
}

/// Nth Catalan number: C(2n, n) / (n+1).
pub fn catalan(n: u64) -> u64 {
    if n == 0 { return 1; }
    combinations(2 * n, n) / (n + 1)
}

/// Bell number B(n) — number of partitions of a set.
pub fn bell(n: u64) -> u64 {
    if n == 0 { return 1; }
    let n = n as usize;
    let mut triangle = vec![vec![0u64; n + 1]; n + 1];
    triangle[0][0] = 1;
    for i in 1..=n {
        triangle[i][0] = triangle[i - 1][i - 1];
        for j in 1..=i {
            triangle[i][j] = triangle[i][j - 1].saturating_add(triangle[i - 1][j - 1]);
        }
    }
    triangle[n][0]
}

/// Stirling number of the second kind S(n, k).
pub fn stirling2(n: u64, k: u64) -> u64 {
    if k == 0 { return if n == 0 { 1 } else { 0 }; }
    if k > n { return 0; }
    if k == 1 || k == n { return 1; }
    let mut result: i64 = 0;
    for i in 0..=k {
        let sign = if (k - i) % 2 == 0 { 1i64 } else { -1i64 };
        let term = combinations(k, i) as i64 * (i as i64).pow(n as u32);
        result += sign * term;
    }
    (result / factorial(k) as i64) as u64
}

/// Stirling number of the first kind |s(n, k)| (unsigned).
pub fn stirling1(n: u64, k: u64) -> u64 {
    if k == 0 { return if n == 0 { 1 } else { 0 }; }
    if k > n { return 0; }
    let n = n as usize;
    let k = k as usize;
    let mut dp = vec![vec![0u64; k + 1]; n + 1];
    dp[0][0] = 1;
    for i in 1..=n {
        for j in 1..=k.min(i) {
            dp[i][j] = dp[i - 1][j - 1] + (i as u64 - 1) * dp[i - 1][j];
        }
    }
    dp[n][k]
}

/// Number of integer partitions of n.
pub fn partition(n: u64) -> u64 {
    let n = n as usize;
    let mut dp = vec![0u64; n + 1];
    dp[0] = 1;
    for i in 1..=n {
        for j in i..=n {
            dp[j] = dp[j].saturating_add(dp[j - i]);
        }
    }
    dp[n]
}

/// Subfactorial (derangements): !n.
pub fn derangements(n: u64) -> u64 {
    if n == 0 { return 1; }
    if n == 1 { return 0; }
    let mut prev2: u64 = 1; // D(0)
    let mut prev1: u64 = 0; // D(1)
    for i in 2..=n {
        let current = (i - 1) * (prev1 + prev2);
        prev2 = prev1;
        prev1 = current;
    }
    prev1
}

/// Nth Fibonacci number.
pub fn fibonacci(n: u64) -> u64 {
    if n == 0 { return 0; }
    if n == 1 { return 1; }
    let mut a: u64 = 0;
    let mut b: u64 = 1;
    for _ in 2..=n {
        let c = a.saturating_add(b);
        a = b;
        b = c;
    }
    b
}

/// Nth Lucas number.
pub fn lucas(n: u64) -> u64 {
    if n == 0 { return 2; }
    if n == 1 { return 1; }
    let mut a: u64 = 2;
    let mut b: u64 = 1;
    for _ in 2..=n {
        let c = a.saturating_add(b);
        a = b;
        b = c;
    }
    b
}

/// Pascal's triangle: returns n rows.
pub fn pascals_triangle(rows: usize) -> Vec<Vec<u64>> {
    let mut triangle: Vec<Vec<u64>> = Vec::with_capacity(rows);
    for i in 0..rows {
        let mut row = vec![1u64; i + 1];
        for j in 1..i {
            row[j] = triangle[i - 1][j - 1] + triangle[i - 1][j];
        }
        triangle.push(row);
    }
    triangle
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_factorial() {
        assert_eq!(factorial(0), 1);
        assert_eq!(factorial(5), 120);
        assert_eq!(factorial(10), 3628800);
    }

    #[test]
    fn test_combinations() {
        assert_eq!(combinations(5, 2), 10);
        assert_eq!(combinations(10, 3), 120);
        assert_eq!(combinations(0, 0), 1);
    }

    #[test]
    fn test_permutations() {
        assert_eq!(permutations(5, 2), 20);
        assert_eq!(permutations(10, 3), 720);
    }

    #[test]
    fn test_fibonacci() {
        assert_eq!(fibonacci(0), 0);
        assert_eq!(fibonacci(1), 1);
        assert_eq!(fibonacci(10), 55);
    }

    #[test]
    fn test_catalan() {
        assert_eq!(catalan(0), 1);
        assert_eq!(catalan(4), 14);
    }

    #[test]
    fn test_bell() {
        assert_eq!(bell(0), 1);
        assert_eq!(bell(1), 1);
        assert_eq!(bell(5), 52);
    }

    #[test]
    fn test_derangements() {
        assert_eq!(derangements(0), 1);
        assert_eq!(derangements(1), 0);
        assert_eq!(derangements(4), 9);
    }

    #[test]
    fn test_partition() {
        assert_eq!(partition(5), 7);
    }

    #[test]
    fn test_lucas() {
        assert_eq!(lucas(0), 2);
        assert_eq!(lucas(1), 1);
        assert_eq!(lucas(5), 11);
    }

    #[test]
    fn test_pascals_triangle() {
        let t = pascals_triangle(5);
        assert_eq!(t[4], vec![1, 4, 6, 4, 1]);
    }
}
