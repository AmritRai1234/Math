<p align="center">
  <h1 align="center">🧮 mathrs</h1>
  <p align="center">
    <strong>A comprehensive, zero-dependency mathematical framework for Rust</strong>
  </p>
  <p align="center">
    <a href="#modules">96 Modules</a> · <a href="#features">10,600+ Lines</a> · <a href="#testing">430+ Tests</a> · <code>no_std</code> Compatible
  </p>
</p>

---

## Overview

**mathrs** is an enterprise-grade, pure-Rust math library that covers virtually every domain of mathematics — from basic arithmetic to abstract algebra, from signal processing to epidemiological modeling. It is designed with three core principles:

- **Zero Dependencies** — No external crates. Everything is implemented from scratch.
- **`no_std` Support** — Works in embedded, WASM, and bare-metal environments (disable the default `std` feature).
- **Comprehensive Testing** — 430+ unit tests ensure correctness across all modules.

## Quick Start

Add `mathrs` to your `Cargo.toml`:

```toml
[dependencies]
mathrs = { git = "https://github.com/AmritRai1234/Math.git" }
```

For `no_std` environments:

```toml
[dependencies]
mathrs = { git = "https://github.com/AmritRai1234/Math.git", default-features = false }
```

## Usage Examples

### Basic Arithmetic & Constants

```rust
use mathrs::constants::*;
use mathrs::arithmetic;

let val = arithmetic::abs(-42.0);        // 42.0
let clamped = arithmetic::clamp(15.0, 0.0, 10.0); // 10.0
println!("π = {}, e = {}, φ = {}", PI, E, PHI);
```

### Linear Algebra

```rust
use mathrs::linalg::*;

let a = Vec3::new(1.0, 2.0, 3.0);
let b = Vec3::new(4.0, 5.0, 6.0);

let dot = a.dot(&b);       // 32.0
let cross = a.cross(&b);   // (-3, 6, -3)
let len = a.magnitude();   // √14
```

### Quaternion Rotations

```rust
use mathrs::quaternion::Quat;

let q = Quat::from_axis_angle(0.0, 1.0, 0.0, std::f64::consts::FRAC_PI_2);
let inv = q.inverse();
let combined = q.mul(&inv); // ≈ identity
```

### Number Theory

```rust
use mathrs::number_theory;

let g = number_theory::gcd(48, 18);          // 6
let prime = number_theory::is_prime(104729);  // true
let tot = number_theory::euler_totient(12);   // 4
```

### Statistics & Distributions

```rust
use mathrs::stats;
use mathrs::distributions;

let data = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
let mean = stats::mean(&data);       // 5.0
let stdev = stats::std_dev(&data);   // ~2.0

let pdf = distributions::normal_pdf(0.0, 0.0, 1.0); // ≈ 0.3989
```

### Signal Processing (FFT)

```rust
use mathrs::signal;

let signal_data = vec![1.0, 0.0, -1.0, 0.0];
// Perform DFT, convolution, windowing, etc.
```

### Calculus

```rust
use mathrs::calculus;

let derivative = calculus::derivative(|x| x * x, 3.0);  // ≈ 6.0
let integral = calculus::integrate(|x| x * x, 0.0, 1.0, 1000); // ≈ 0.333
```

### Graph Theory

```rust
use mathrs::graph::Graph;

let mut g = Graph::new(5);
g.add_edge(0, 1, 1.0);
g.add_edge(1, 2, 2.0);
let dist = g.dijkstra(0); // Shortest paths from node 0
```

## Modules

**mathrs** ships with **96 modules** organized across every major branch of mathematics:

### 🔢 Core Mathematics
| Module | Description |
|--------|-------------|
| `constants` | Mathematical constants (π, e, φ, √2, etc.) |
| `arithmetic` | Absolute value, rounding, clamping, lerp, remapping |
| `powers` | Powers, roots, sqrt, cbrt, hypot |
| `exponential` | exp, ln, log, sigmoid, softplus, logit |
| `trig` | sin, cos, tan, asin, acos, atan, atan2 |
| `hyperbolic` | sinh, cosh, tanh, asinh, acosh, atanh |
| `number_theory` | GCD, LCM, primes, modular arithmetic, Euler's totient |
| `number_bases` | Base conversion (binary, octal, hex, arbitrary) |
| `sequences` | Fibonacci, Lucas, Padovan, tribonacci, A/GP series |
| `combinatorics` | Factorial, permutations, combinations, Catalan, Bell |
| `polynomial` | Polynomial arithmetic, evaluation, roots, derivatives |

### 📐 Algebra & Structures
| Module | Description |
|--------|-------------|
| `abstract_algebra` | Groups, rings, fields, homomorphisms |
| `linalg` | Vec2/3/4, Mat2/3/4, transforms, projections |
| `matrix` | General MxN matrix operations, determinant, inverse |
| `quaternion` | Quaternion algebra, rotations, SLERP |
| `complex` | Complex number arithmetic, polar form |
| `rational` | Rational number arithmetic (exact fractions) |
| `bigint` | Arbitrary-precision integer arithmetic |
| `tensor` | Multi-dimensional tensor operations |
| `sparse` | Sparse matrix storage and operations |
| `finite_fields` | Galois field arithmetic (GF(p)) |
| `tropical` | Tropical (min-plus / max-plus) algebra |

### 📊 Statistics & Probability
| Module | Description |
|--------|-------------|
| `stats` | Mean, median, variance, correlation, percentiles |
| `distributions` | Normal, binomial, Poisson, chi-squared, exponential |
| `probability` | Conditional probability, Bayes' theorem, entropy |
| `bayesian` | Bayesian inference and updating |
| `hypothesis` | Hypothesis testing (z-test, t-test, chi-squared) |
| `regression` | Linear, polynomial, and logistic regression |
| `stochastic` | Stochastic processes, Brownian motion |
| `markov` | Markov chains, transition matrices, steady state |
| `random` | Pseudorandom number generation (LCG, xorshift) |

### 📈 Analysis & Calculus
| Module | Description |
|--------|-------------|
| `calculus` | Numerical derivatives, integrals, gradient, Jacobian |
| `rootfinding` | Bisection, Newton-Raphson, Brent's method |
| `ode` | ODE solvers: Euler, RK4, RK45, Verlet |
| `pde` | Partial differential equation solvers |
| `numerical` | General numerical methods and solvers |
| `autodiff` | Automatic differentiation (forward mode) |
| `optimization` | Gradient descent, Nelder-Mead, golden section |
| `equations` | Equation solving and systems of equations |
| `inequalities` | Inequality verification (AM-GM, Cauchy-Schwarz) |
| `error_analysis` | Error propagation and numerical error bounds |
| `continued_fractions` | Continued fraction expansions and convergents |
| `dynamical_systems` | Lorenz, logistic map, Hénon, phase portraits |

### 🌐 Geometry & Topology
| Module | Description |
|--------|-------------|
| `geometry` | Area, volume, distances, triangles, polygons |
| `coordinates` | Polar, cylindrical, spherical coordinate transforms |
| `coordinate` | Additional coordinate system utilities |
| `curves` | Parametric curves, arc length, curvature |
| `spline` | Cubic spline interpolation |
| `convex` | Convex hull, containment, support functions |
| `diff_geometry` | Differential geometry: curvature, torsion, geodesics |
| `topology` | Euler characteristic, simplicial complexes, Betti numbers |
| `geospatial` | Haversine, Vincenty, map projections |
| `fractals` | Mandelbrot, Julia, Sierpinski, Koch, IFS |

### 📡 Signal & Data Processing
| Module | Description |
|--------|-------------|
| `signal` | FFT, DFT, DCT, convolution, windowing |
| `signal_processing` | Filters, spectral analysis, resampling |
| `wavelets` | Haar, Daubechies wavelet transforms |
| `noise` | Perlin, simplex, fBm noise generation |
| `interpolation` | Lerp, Bézier, Hermite, Lagrange, smoothstep |
| `time_series` | Moving averages, exponential smoothing, trends |
| `compression` | Run-length encoding, Huffman, LZ compression |
| `image_processing` | Convolution kernels, edge detection, filters |
| `clustering` | K-means, DBSCAN, hierarchical clustering |
| `dimensionality` | PCA, dimensionality reduction techniques |

### 🔐 Information & Coding Theory
| Module | Description |
|--------|-------------|
| `information` | Shannon entropy, mutual information |
| `information_theory` | Channel capacity, source coding theorems |
| `coding` | Error-correcting codes: Hamming, Reed-Solomon |
| `crypto_math` | Modular exponentiation, RSA primitives, hashing |
| `string_metrics` | Levenshtein, Hamming distance, Jaccard similarity |

### 🎮 Applied Mathematics
| Module | Description |
|--------|-------------|
| `physics` | Kinematics, gravitation, optics, thermodynamics |
| `financial` | Compound interest, NPV, Black-Scholes, amortization |
| `game_theory` | Nash equilibrium, minimax, payoff matrices |
| `voting` | Voting systems (plurality, Borda, Condorcet) |
| `epidemiology` | SIR/SIS/SEIR compartmental models |
| `control_theory` | PID controllers, transfer functions, stability |
| `queueing` | M/M/1, M/M/c queueing models, Little's law |
| `date_math` | Julian day, day-of-week, leap year calculations |
| `color_math` | RGB/HSL/HSV conversion, color blending, contrast |
| `music_theory` | Frequency ratios, intervals, MIDI conversion |
| `recreational` | Magic squares, Tower of Hanoi, happy numbers |

### 🧠 Logic & Computation
| Module | Description |
|--------|-------------|
| `logic` | Boolean algebra, truth tables, logical operators |
| `logic_prop` | Propositional logic, CNF, satisfiability |
| `sets` | Set operations: union, intersection, power set |
| `set_theory` | Formal set theory operations and relations |
| `bits` | Bit manipulation: popcount, clz, rotate, parity |
| `automata` | DFA, NFA, regex engines, state machines |
| `cellular_automata` | Elementary CA, Game of Life, rule systems |
| `fuzzy` | Fuzzy logic: membership, operations, defuzzification |
| `order_theory` | Partial orders, lattices, Hasse diagrams |

### 🔬 Advanced Topics
| Module | Description |
|--------|-------------|
| `special` | Gamma, beta, Bessel, zeta, erf, Lambert W, Airy |
| `interval` | Interval arithmetic with guaranteed bounds |
| `decomposition` | LU, QR, Cholesky, SVD matrix decomposition |
| `graph` | Graph theory: BFS, DFS, Dijkstra, MST |
| `network_flow` | Max-flow, min-cut, Ford-Fulkerson |
| `diophantine` | Diophantine equation solvers |
| `units` | Unit conversion (length, mass, temperature, etc.) |

## Project Structure

```
Math/
├── Cargo.toml          # Package manifest (zero dependencies)
├── Cargo.lock
├── .gitignore
├── README.md
└── src/
    ├── lib.rs           # Module declarations & crate docs
    ├── constants.rs     # Mathematical constants
    ├── arithmetic.rs    # Basic arithmetic operations
    ├── linalg.rs        # Linear algebra (vectors & matrices)
    ├── ... (96 modules)
    └── decomposition.rs # Matrix decompositions
```

## Testing

Run the full test suite:

```bash
cargo test
```

Run tests for a specific module:

```bash
cargo test --lib quaternion
cargo test --lib number_theory
cargo test --lib stats
```

All **430+ tests** cover correctness with floating-point tolerance where applicable.

## Features

| Feature | Default | Description |
|---------|---------|-------------|
| `std` | ✅ | Enables standard library support |

Disable `std` for `no_std` / embedded / WASM targets:

```toml
[dependencies]
mathrs = { git = "https://github.com/AmritRai1234/Math.git", default-features = false }
```

## Design Philosophy

1. **From-Scratch Implementations** — Every algorithm is hand-written in Rust. No hidden dependencies, no FFI, no C bindings.
2. **Inline & Const-Friendly** — Hot-path functions are `#[inline]` annotated for optimal codegen.
3. **Modular Architecture** — Each mathematical domain lives in its own module. Import only what you need.
4. **Correctness First** — Numerical accuracy is validated against known mathematical identities and edge cases.

## Requirements

- **Rust** 1.56+ (2021 edition)
- No external dependencies

## License

This project is licensed under the **MIT License**.

---

<p align="center">
  Built with ❤️ and pure math
</p>
