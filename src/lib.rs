//! # mathrs
//!
//! A comprehensive, zero-dependency math library for Rust.
//!
//! ## Modules
//!
//! - **constants** — Mathematical constants (π, e, φ, etc.)
//! - **arithmetic** — Basic arithmetic, rounding, clamping, lerp
//! - **powers** — Powers, roots, sqrt, cbrt, hypot
//! - **exponential** — exp, ln, log, sigmoid, softplus
//! - **trig** — Trigonometric functions (sin, cos, tan, asin, etc.)
//! - **hyperbolic** — Hyperbolic functions (sinh, cosh, tanh, etc.)
//! - **number_theory** — GCD, primes, modular arithmetic, totient
//! - **combinatorics** — Factorial, combinations, Fibonacci, Catalan, Bell
//! - **stats** — Mean, median, variance, correlation, percentiles
//! - **distributions** — Normal, binomial, Poisson, chi-squared, etc.
//! - **special** — Gamma, beta, Bessel, zeta, erf, Lambert W, Airy
//! - **linalg** — Vectors (Vec2/3/4), matrices (Mat2/3/4), transforms
//! - **quaternion** — Quaternion math, rotations, slerp
//! - **complex** — Complex number arithmetic and functions
//! - **interpolation** — Lerp, Bézier, Hermite, Lagrange, smoothstep
//! - **calculus** — Numerical derivatives, integrals, gradient, Jacobian
//! - **rootfinding** — Bisection, Newton-Raphson, Brent, golden section
//! - **ode** — ODE solvers: Euler, RK4, RK45, Verlet
//! - **signal** — FFT, DFT, DCT, convolution, windowing
//! - **noise** — Perlin, simplex, fBm, RNG
//! - **sets** — Set operations: union, intersection, power set
//! - **bits** — Bit manipulation: popcount, clz, rotate, parity

pub mod constants;
pub mod arithmetic;
pub mod exponential;
pub mod powers;
pub mod trig;
pub mod hyperbolic;
pub mod number_theory;
pub mod combinatorics;
pub mod stats;
pub mod distributions;
pub mod special;
pub mod linalg;
pub mod quaternion;
pub mod complex;
pub mod interpolation;
pub mod calculus;
pub mod rootfinding;
pub mod ode;
pub mod signal;
pub mod noise;
pub mod sets;
pub mod bits;
pub mod polynomial;
pub mod coordinates;
pub mod geometry;
pub mod sequences;
pub mod financial;
pub mod matrix;
pub mod rational;
pub mod graph;
pub mod regression;
pub mod information;
pub mod continued_fractions;
pub mod logic;
pub mod crypto_math;
pub mod game_theory;
pub mod markov;
pub mod tensor;
pub mod pde;
pub mod units;
pub mod bigint;
pub mod curves;
pub mod optimization;
pub mod sparse;
pub mod clustering;
pub mod time_series;
pub mod wavelets;
pub mod autodiff;
pub mod interval;
pub mod fuzzy;
pub mod coding;
pub mod abstract_algebra;
pub mod queueing;
pub mod voting;
pub mod random;
pub mod spline;
pub mod fractals;
pub mod stochastic;
pub mod automata;
pub mod topology;
pub mod physics;
pub mod number_bases;
pub mod tropical;
pub mod error_analysis;
pub mod cellular_automata;
pub mod string_metrics;
pub mod dimensionality;
pub mod diff_geometry;
pub mod geospatial;
pub mod color_math;
pub mod music_theory;
pub mod hypothesis;
pub mod control_theory;
pub mod compression;
pub mod network_flow;
pub mod date_math;
pub mod image_processing;
pub mod probability;
pub mod bayesian;
pub mod order_theory;
pub mod signal_processing;
pub mod coordinate;
pub mod set_theory;
pub mod equations;
pub mod recreational;
pub mod epidemiology;
pub mod inequalities;
pub mod finite_fields;
pub mod numerical;
pub mod information_theory;
pub mod convex;
pub mod dynamical_systems;
pub mod diophantine;
pub mod decomposition;
pub mod logic_prop;
