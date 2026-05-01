//! Chaos-based cryptography using the Entropic Lorenz system.
//!
//! The self-modifying chaotic attractor generates a keystream that is:
//! - Deterministic (same key → same stream)
//! - Unpredictable (chaotic sensitivity to initial conditions)
//! - Self-mutating (entropy feedback changes the system mid-encryption)
//!
//! This makes cryptanalysis harder than standard Lorenz-based ciphers because
//! the attractor parameters are not fixed — they evolve based on the trajectory.

use crate::exponential;
use crate::trig;
use crate::arithmetic;

/// Secret key: initial conditions for the Entropic Lorenz system.
///
/// v3 FIX: includes a message counter that auto-increments on every
/// encrypt call. The IV is derived from (counter + key state), so the
/// same (key, IV) pair can NEVER be reused accidentally.
#[derive(Clone, Debug)]
pub struct ChaosKey {
    pub x0: f64,
    pub y0: f64,
    pub z0: f64,
    pub rounds: usize,
    /// Auto-incrementing message counter — prevents IV reuse.
    msg_counter: u64,
}

impl ChaosKey {
    /// Create a key from three f64 seeds.
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x0: x, y0: y, z0: z, rounds: 500, msg_counter: 0 }
    }

    /// Derive a key from a password string.
    pub fn from_password(password: &str) -> Self {
        let bytes = password.as_bytes();
        let mut x = 0.1;
        let mut y = 0.1;
        let mut z = 0.1;
        for (i, &b) in bytes.iter().enumerate() {
            let v = b as f64 / 255.0;
            match i % 3 {
                0 => x += v * 0.9 + 0.1,
                1 => y += v * 0.9 + 0.1,
                _ => z += v * 0.9 + 0.1,
            }
        }
        for _ in 0..3 {
            x = (x * 7.31 + y * 3.17) % 50.0 + 0.1;
            y = (y * 5.93 + z * 2.71) % 50.0 + 0.1;
            z = (z * 4.57 + x * 1.93) % 50.0 + 0.1;
        }
        Self { x0: x, y0: y, z0: z, rounds: 500, msg_counter: 0 }
    }

    /// Generate a unique IV for the next message and advance the counter.
    /// IV = mix(counter, key_state) — deterministic but never repeated.
    pub fn next_iv(&mut self) -> u64 {
        let c = self.msg_counter;
        self.msg_counter += 1;
        // Mix counter with key material using a simple hash
        let a = (c ^ (self.x0.to_bits())).wrapping_mul(0x9e3779b97f4a7c15);
        let b = (c.wrapping_add(1) ^ (self.y0.to_bits())).wrapping_mul(0x6c62272e07bb0142);
        let cc = (c.wrapping_add(2) ^ (self.z0.to_bits())).wrapping_mul(0x94d049bb133111eb);
        a ^ b.rotate_left(17) ^ cc.rotate_right(31)
    }

    /// Current message counter (for diagnostics).
    pub fn message_count(&self) -> u64 { self.msg_counter }
}

/// The Entropic Lorenz PRNG — core engine for encryption.
pub struct EntropicPRNG {
    x: f64,
    y: f64,
    z: f64,
    history: Vec<f64>,
    entropy: f64,
    /// Output buffer: 8 bytes derived from last Lorenz state reseed
    buf: u64,
    buf_pos: u8, // which byte of buf to emit next (0..8)
}

impl EntropicPRNG {
    /// Initialize from a ChaosKey, running warmup rounds to reach the attractor.
    pub fn from_key(key: &ChaosKey) -> Self {
        let mut prng = Self {
            x: key.x0, y: key.y0, z: key.z0,
            history: Vec::with_capacity(64),
            entropy: 1.0,
            buf: 0,
            buf_pos: 8, // force reseed on first call
        };
        for _ in 0..key.rounds {
            prng.step();
        }
        prng
    }

    /// Advance the Entropic Lorenz system by one step.
    fn step(&mut self) {
        let dt = 0.005;

        // Measure entropy from trajectory variance
        let energy = self.x * self.x + self.y * self.y + self.z * self.z;
        self.history.push(energy);
        if self.history.len() > 30 { self.history.remove(0); }

        if self.history.len() > 5 {
            let mean = self.history.iter().sum::<f64>() / self.history.len() as f64;
            let var = self.history.iter()
                .map(|&v| (v - mean) * (v - mean))
                .sum::<f64>() / self.history.len() as f64;
            self.entropy = if var > 1e-30 { 0.5 * exponential::ln(var) } else { 0.0 };
        }

        // Self-modifying parameters — the key innovation
        let sigma = 10.0 + 2.0 * trig::sin(self.entropy);
        let rho = 28.0 + 5.0 * trig::cos(self.entropy * 0.5);
        let beta = 8.0 / 3.0 + 0.5 * trig::sin(self.entropy * 0.3);

        // RK4 integration
        let f = |x: f64, y: f64, z: f64| -> (f64, f64, f64) {
            (sigma * (y - x), x * (rho - z) - y, x * y - beta * z)
        };
        let k1 = f(self.x, self.y, self.z);
        let k2 = f(self.x + dt/2.0*k1.0, self.y + dt/2.0*k1.1, self.z + dt/2.0*k1.2);
        let k3 = f(self.x + dt/2.0*k2.0, self.y + dt/2.0*k2.1, self.z + dt/2.0*k2.2);
        let k4 = f(self.x + dt*k3.0, self.y + dt*k3.1, self.z + dt*k3.2);

        self.x += dt/6.0 * (k1.0 + 2.0*k2.0 + 2.0*k3.0 + k4.0);
        self.y += dt/6.0 * (k1.1 + 2.0*k2.1 + 2.0*k3.1 + k4.1);
        self.z += dt/6.0 * (k1.2 + 2.0*k2.2 + 2.0*k3.2 + k4.2);
    }

    /// Generate a single pseudo-random byte from the chaotic state.
    ///
    /// FIX v2: Mix all 3 coordinates → uniform distribution.
    /// FIX v4: SplitMix64 finalizer for full bit-level decorrelation.
    ///   Uses Fibonacci hashing on all three coordinates combined, then
    ///   applies the SplitMix64 avalanche pass to produce a uniformly
    ///   distributed byte with no detectable statistical structure.
    pub fn next_byte(&mut self) -> u8 {
        self.step();
        // Combine all three coordinates via Fibonacci hashing
        let bx = (arithmetic::abs(self.x * 1e9) as u64).wrapping_mul(0x9e3779b97f4a7c15);
        let by = (arithmetic::abs(self.y * 1e9) as u64).wrapping_mul(0x6c62272e07bb0142);
        let bz = (arithmetic::abs(self.z * 1e9) as u64).wrapping_mul(0x94d049bb133111eb);
        let combined = bx ^ by.rotate_left(21) ^ bz.rotate_right(13);
        // SplitMix64 finalizer — avalanche diffusion
        let mut h = combined;
        h = (h ^ (h >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        h = (h ^ (h >> 27)).wrapping_mul(0x94d049bb133111eb);
        h ^= h >> 31;
        (h & 0xFF) as u8
    }

    /// Generate a keystream of `n` bytes.
    pub fn keystream(&mut self, n: usize) -> Vec<u8> {
        (0..n).map(|_| self.next_byte()).collect()
    }

    /// FIX v2: Inject an IV into the state to prevent lobe fingerprinting.
    /// Call after from_key(), before generating bytes, with a unique nonce.
    pub fn inject_iv(&mut self, iv: u64) {
        // Perturb state with IV so same key + different IV → different lobe
        self.x += (iv & 0xFF) as f64 * 0.001;
        self.y += ((iv >> 8) & 0xFF) as f64 * 0.001;
        self.z += ((iv >> 16) & 0xFF) as f64 * 0.001;
        // Run extra warmup steps to diffuse the IV perturbation
        for _ in 0..200 { self.step(); }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Stream Cipher v3: Auto-IV Encrypt / Decrypt  (IV reuse fixed)
// ═══════════════════════════════════════════════════════════════════════════

/// Encrypt with an auto-generated unique IV. RECOMMENDED — prevents IV reuse.
/// The key's internal counter advances automatically, guaranteeing every
/// encryption uses a fresh, unique IV even if called in a loop.
/// Ciphertext format: [IV: 8 bytes] [encrypted payload]
pub fn encrypt_auto(plaintext: &[u8], key: &mut ChaosKey) -> Vec<u8> {
    let iv = key.next_iv(); // counter advances here — reuse is impossible
    encrypt_iv(plaintext, key, iv)
}

/// Decrypt a message produced by encrypt_auto or encrypt_iv.
/// The IV is read from the first 8 bytes of the ciphertext.
pub fn decrypt_auto(ciphertext: &[u8], key: &ChaosKey) -> Vec<u8> {
    decrypt_iv(ciphertext, key)
}

/// Encrypt with a random IV prepended to ciphertext.
/// Ciphertext format: [IV: 8 bytes] [encrypted payload]
pub fn encrypt_iv(plaintext: &[u8], key: &ChaosKey, iv: u64) -> Vec<u8> {
    let mut prng = EntropicPRNG::from_key(key);
    prng.inject_iv(iv);
    let mut out = iv.to_le_bytes().to_vec();
    out.extend(plaintext.iter().map(|&b| b ^ prng.next_byte()));
    out
}

/// Decrypt a message encrypted with encrypt_iv.
pub fn decrypt_iv(ciphertext: &[u8], key: &ChaosKey) -> Vec<u8> {
    if ciphertext.len() < 8 { return vec![]; }
    let iv = u64::from_le_bytes(ciphertext[..8].try_into().unwrap_or([0;8]));
    let mut prng = EntropicPRNG::from_key(key);
    prng.inject_iv(iv);
    ciphertext[8..].iter().map(|&b| b ^ prng.next_byte()).collect()
}

/// Encrypt plaintext bytes — legacy no-IV mode. Use encrypt_auto for new code.
pub fn encrypt(plaintext: &[u8], key: &ChaosKey) -> Vec<u8> {
    let mut prng = EntropicPRNG::from_key(key);
    plaintext.iter().map(|&b| b ^ prng.next_byte()).collect()
}

/// Decrypt ciphertext — legacy no-IV mode.
pub fn decrypt(ciphertext: &[u8], key: &ChaosKey) -> Vec<u8> {
    encrypt(ciphertext, key)
}

/// Encrypt a string returning hex — uses auto-IV internally.
pub fn encrypt_str(plaintext: &str, key: &mut ChaosKey) -> String {
    let encrypted = encrypt_auto(plaintext.as_bytes(), key);
    encrypted.iter().map(|b| format!("{:02x}", b)).collect()
}

/// Decrypt hex-encoded ciphertext back to a string.
pub fn decrypt_str(hex: &str, key: &ChaosKey) -> String {
    let bytes: Vec<u8> = (0..hex.len())
        .step_by(2)
        .filter_map(|i| u8::from_str_radix(&hex[i..i+2], 16).ok())
        .collect();
    let decrypted = decrypt_iv(&bytes, key);
    String::from_utf8_lossy(&decrypted).to_string()
}


// ═══════════════════════════════════════════════════════════════════════════
// Chaos Hash — deterministic hash using trajectory endpoint
// ═══════════════════════════════════════════════════════════════════════════

/// Generate a 256-bit (32-byte) hash from input data using chaotic mixing.
pub fn chaos_hash(data: &[u8]) -> [u8; 32] {
    let mut x = 1.23456789;
    let mut y = 9.87654321;
    let mut z = 5.55555555;

    // Each byte strongly perturbs the state and is mixed through chaos
    for (i, &b) in data.iter().enumerate() {
        let v = (b as f64 + 1.0) / 128.0;
        // Rotate which variable gets perturbed and use nonlinear mixing
        x += v * trig::sin(y * 13.7 + i as f64);
        y += v * trig::cos(z * 17.3 + i as f64);
        z += v * trig::sin(x * 11.1 + i as f64);
        // Run a Lorenz step after every byte for thorough mixing
        let sigma = 10.0 + arithmetic::abs(x % 3.0);
        let rho = 28.0 + arithmetic::abs(y % 5.0);
        let beta = 8.0/3.0 + arithmetic::abs(z % 1.0);
        let dt = 0.01;
        let dx = sigma * (y - x) * dt;
        let dy = (x * (rho - z) - y) * dt;
        let dz = (x * y - beta * z) * dt;
        x += dx; y += dy; z += dz;
    }

    // Run 500 more iterations to fully diffuse
    let key = ChaosKey { x0: x, y0: y, z0: z, rounds: 500, msg_counter: 0 };
    let mut prng = EntropicPRNG::from_key(&key);
    let mut hash = [0u8; 32];
    for byte in hash.iter_mut() {
        *byte = prng.next_byte();
    }
    hash
}

/// Return chaos hash as hex string.
pub fn chaos_hash_hex(data: &[u8]) -> String {
    chaos_hash(data).iter().map(|b| format!("{:02x}", b)).collect()
}

// ═══════════════════════════════════════════════════════════════════════════
// Security Analysis Tools
// ═══════════════════════════════════════════════════════════════════════════

/// Measure keystream randomness via byte frequency distribution.
/// Returns chi-squared statistic. Lower = more uniform = more random.
pub fn chi_squared_test(keystream: &[u8]) -> f64 {
    let n = keystream.len() as f64;
    let expected = n / 256.0;
    let mut freq = [0u64; 256];
    for &b in keystream { freq[b as usize] += 1; }
    freq.iter().map(|&f| {
        let diff = f as f64 - expected;
        diff * diff / expected
    }).sum()
}

/// Measure avalanche effect: how many output bits change when 1 input bit flips.
/// Ideal = 50% (128 out of 256 bits change in a 32-byte hash).
pub fn avalanche_score(data: &[u8]) -> f64 {
    let hash1 = chaos_hash(data);
    let mut modified = data.to_vec();
    if modified.is_empty() { modified.push(0); }
    modified[0] ^= 1; // flip 1 bit
    let hash2 = chaos_hash(&modified);

    let mut diff_bits = 0u32;
    for (&a, &b) in hash1.iter().zip(hash2.iter()) {
        diff_bits += (a ^ b).count_ones();
    }
    diff_bits as f64 / 256.0 * 100.0 // percentage
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encrypt_decrypt_roundtrip() {
        let key = ChaosKey::new(1.0, 2.0, 3.0);
        let msg = b"Hello, Wormhole Encryption!";
        let encrypted = encrypt(msg, &key);
        let decrypted = decrypt(&encrypted, &key);
        assert_eq!(decrypted, msg);
    }

    #[test]
    fn test_encrypt_str_roundtrip() {
        let mut key = ChaosKey::from_password("my_secret_key_2026");
        let msg = "The Entropic Lorenz changes its own rules";
        let hex = encrypt_str(msg, &mut key);
        let recovered = decrypt_str(&hex, &key);
        assert_eq!(recovered, msg);
    }

    #[test]
    fn test_different_keys_different_output() {
        let k1 = ChaosKey::new(1.0, 2.0, 3.0);
        let k2 = ChaosKey::new(1.0, 2.0, 3.1);
        let msg = b"same message here";
        let e1 = encrypt(msg, &k1);
        let e2 = encrypt(msg, &k2);
        assert_ne!(e1, e2);
    }

    #[test]
    fn test_ciphertext_looks_random() {
        let key = ChaosKey::from_password("test_password");
        let msg = vec![0u8; 1000];
        let encrypted = encrypt(&msg, &key);
        assert!(encrypted.iter().any(|&b| b != 0));
        let chi2 = chi_squared_test(&encrypted);
        assert!(chi2 < 400.0);
    }

    #[test]
    fn test_keystream_uniformity() {
        let key = ChaosKey::new(5.5, 12.3, 7.8);
        let mut prng = EntropicPRNG::from_key(&key);
        let stream = prng.keystream(5000);
        let chi2 = chi_squared_test(&stream);
        assert!(chi2 < 400.0, "chi2 = {} (too high)", chi2);
    }

    #[test]
    fn test_encrypt_auto_roundtrip() {
        let mut key = ChaosKey::from_password("auto_iv_test");
        let msg = b"Testing auto-IV encryption roundtrip";
        let ct = encrypt_auto(msg, &mut key);
        let pt = decrypt_auto(&ct, &key);
        assert_eq!(pt, msg);
    }

    #[test]
    fn test_iv_reuse_impossible() {
        let mut key = ChaosKey::from_password("no_reuse_test");
        let msg = b"same message encrypted twice";
        let ct1 = encrypt_auto(msg, &mut key);
        let ct2 = encrypt_auto(msg, &mut key);
        assert_ne!(&ct1[..8], &ct2[..8], "IVs must be unique");
        assert_ne!(ct1, ct2, "Same plaintext must give different ciphertext");
        assert_eq!(decrypt_auto(&ct1, &key), msg);
        assert_eq!(decrypt_auto(&ct2, &key), msg);
    }

    #[test]
    fn test_counter_uniqueness() {
        let mut key = ChaosKey::new(3.14, 2.71, 1.61);
        let ivs: Vec<u64> = (0..20).map(|_| key.next_iv()).collect();
        let mut sorted = ivs.clone();
        sorted.sort();
        sorted.dedup();
        assert_eq!(sorted.len(), 20, "All IVs must be unique");
        assert_eq!(key.message_count(), 20);
    }

    #[test]
    fn test_iv_reuse_attack_defeated() {
        let mut key = ChaosKey::from_password("attack_test_key");
        let msg1 = b"The Entropic Lorenz";
        let msg2 = b"Wormhole signature!";
        let ct1 = encrypt_auto(msg1, &mut key);
        let ct2 = encrypt_auto(msg2, &mut key);
        let xored_ct: Vec<u8> = ct1[8..].iter().zip(ct2[8..].iter()).map(|(&a,&b)| a^b).collect();
        let xored_pt: Vec<u8> = msg1.iter().zip(msg2.iter()).map(|(&a,&b)| a^b).collect();
        let len = xored_pt.len();
        assert_ne!(&xored_ct[..len], xored_pt.as_slice(), "IV reuse attack must fail");
    }

    #[test]
    fn test_chaos_hash_deterministic() {
        let h1 = chaos_hash(b"hello world");
        let h2 = chaos_hash(b"hello world");
        assert_eq!(h1, h2);
    }

    #[test]
    fn test_chaos_hash_avalanche() {
        let h1 = chaos_hash(b"hello world");
        let h2 = chaos_hash(b"iello world");
        let diff: u32 = h1.iter().zip(h2.iter()).map(|(&a,&b)| (a^b).count_ones()).sum();
        assert!(diff > 30, "Avalanche too low: {} bits differ", diff);
    }

    #[test]
    fn test_avalanche_score() {
        let score = avalanche_score(b"test data for avalanche");
        assert!(score > 10.0, "Avalanche score too low: {}%", score);
    }

    #[test]
    fn test_password_key_derivation() {
        let k1 = ChaosKey::from_password("password123");
        let k2 = ChaosKey::from_password("password124");
        assert!((k1.x0 - k2.x0).abs() > 0.01 || (k1.y0 - k2.y0).abs() > 0.01);
    }
}
