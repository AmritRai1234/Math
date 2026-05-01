//! Physics formulas: kinematics, electromagnetism, thermodynamics, orbital mechanics.

use crate::powers;
use crate::trig;
use crate::constants::PI;
use crate::units;

// ─── Kinematics ───

/// Position under constant acceleration: x = x0 + v0*t + 0.5*a*t².
pub fn position(x0:f64,v0:f64,a:f64,t:f64)->f64{x0+v0*t+0.5*a*t*t}
/// Velocity under constant acceleration: v = v0 + a*t.
pub fn velocity(v0:f64,a:f64,t:f64)->f64{v0+a*t}
/// v² = v0² + 2*a*(x - x0).
pub fn velocity_from_displacement(v0:f64,a:f64,dx:f64)->f64{powers::sqrt(v0*v0+2.0*a*dx)}
/// Time of flight for projectile.
pub fn projectile_time(v0:f64,angle:f64,g:f64)->f64{2.0*v0*trig::sin(angle)/g}
/// Max height of projectile.
pub fn projectile_max_height(v0:f64,angle:f64,g:f64)->f64{let vs=v0*trig::sin(angle);vs*vs/(2.0*g)}
/// Range of projectile.
pub fn projectile_range(v0:f64,angle:f64,g:f64)->f64{v0*v0*trig::sin(2.0*angle)/g}

// ─── Newton's Laws ───

/// Force = mass * acceleration.
pub fn force(m:f64,a:f64)->f64{m*a}
/// Gravitational force between two masses.
pub fn gravitational_force(m1:f64,m2:f64,r:f64)->f64{units::GRAVITATIONAL*m1*m2/(r*r)}
/// Gravitational potential energy.
pub fn gravitational_pe(m1:f64,m2:f64,r:f64)->f64{-units::GRAVITATIONAL*m1*m2/r}
/// Kinetic energy.
pub fn kinetic_energy(m:f64,v:f64)->f64{0.5*m*v*v}
/// Potential energy (near surface).
pub fn potential_energy(m:f64,g:f64,h:f64)->f64{m*g*h}
/// Work = F * d * cos(θ).
pub fn work(f:f64,d:f64,theta:f64)->f64{f*d*trig::cos(theta)}
/// Power = work / time.
pub fn power(w:f64,t:f64)->f64{w/t}
/// Momentum.
pub fn momentum(m:f64,v:f64)->f64{m*v}
/// Impulse = F * Δt.
pub fn impulse(f:f64,dt:f64)->f64{f*dt}

// ─── Rotational ───

/// Centripetal acceleration: a = v²/r.
pub fn centripetal_accel(v:f64,r:f64)->f64{v*v/r}
/// Angular velocity: ω = v/r.
pub fn angular_velocity(v:f64,r:f64)->f64{v/r}
/// Moment of inertia (point mass).
pub fn moment_of_inertia_point(m:f64,r:f64)->f64{m*r*r}
/// Torque = r × F.
pub fn torque(r:f64,f:f64,theta:f64)->f64{r*f*trig::sin(theta)}

// ─── Electromagnetism ───

/// Coulomb's law: F = k*q1*q2/r².
pub fn coulomb_force(q1:f64,q2:f64,r:f64)->f64{8.9875e9*q1*q2/(r*r)}
/// Electric field from point charge.
pub fn electric_field(q:f64,r:f64)->f64{8.9875e9*q/(r*r)}
/// Electric potential.
pub fn electric_potential(q:f64,r:f64)->f64{8.9875e9*q/r}
/// Capacitance of parallel plates.
pub fn capacitance_parallel(area:f64,d:f64,epsilon:f64)->f64{epsilon*area/d}
/// Ohm's law: V = I*R.
pub fn ohms_law_v(i:f64,r:f64)->f64{i*r}
pub fn ohms_law_i(v:f64,r:f64)->f64{v/r}
pub fn ohms_law_r(v:f64,i:f64)->f64{v/i}
/// Resistors in series.
pub fn resistance_series(resistors:&[f64])->f64{resistors.iter().sum()}
/// Resistors in parallel.
pub fn resistance_parallel(resistors:&[f64])->f64{1.0/resistors.iter().map(|&r|1.0/r).sum::<f64>()}
/// Magnetic force on moving charge: F = qvB sin(θ).
pub fn lorentz_force(q:f64,v:f64,b:f64,theta:f64)->f64{q*v*b*trig::sin(theta)}

// ─── Waves ───

/// Wave speed: v = f*λ.
pub fn wave_speed(freq:f64,wavelength:f64)->f64{freq*wavelength}
/// Frequency from period.
pub fn frequency(period:f64)->f64{1.0/period}
/// Doppler effect (source moving toward observer).
pub fn doppler(f0:f64,v_sound:f64,v_source:f64)->f64{f0*v_sound/(v_sound-v_source)}
/// De Broglie wavelength: λ = h/(mv).
pub fn de_broglie(m:f64,v:f64)->f64{units::PLANCK/(m*v)}

// ─── Thermodynamics ───

/// Ideal gas law: PV = nRT. Solve for P.
pub fn ideal_gas_pressure(n:f64,t:f64,v:f64)->f64{n*units::GAS_CONSTANT*t/v}
/// Entropy change: ΔS = Q/T.
pub fn entropy_change(q:f64,t:f64)->f64{q/t}
/// Carnot efficiency.
pub fn carnot_efficiency(t_hot:f64,t_cold:f64)->f64{1.0-t_cold/t_hot}
/// Stefan-Boltzmann radiation: P = σ*A*T⁴.
pub fn stefan_boltzmann(area:f64,temp:f64)->f64{5.670374419e-8*area*powers::pow(temp,4.0)}

// ─── Orbital Mechanics ───

/// Orbital velocity: v = sqrt(GM/r).
pub fn orbital_velocity(m_central:f64,r:f64)->f64{powers::sqrt(units::GRAVITATIONAL*m_central/r)}
/// Escape velocity: v = sqrt(2GM/r).
pub fn escape_velocity(m:f64,r:f64)->f64{powers::sqrt(2.0*units::GRAVITATIONAL*m/r)}
/// Orbital period: T = 2π*sqrt(r³/(GM)).
pub fn orbital_period(m:f64,r:f64)->f64{2.0*PI*powers::sqrt(r*r*r/(units::GRAVITATIONAL*m))}
/// Schwarzschild radius: r_s = 2GM/c².
pub fn schwarzschild_radius(m:f64)->f64{2.0*units::GRAVITATIONAL*m/(units::SPEED_OF_LIGHT*units::SPEED_OF_LIGHT)}

// ─── Relativity ───

/// Lorentz factor: γ = 1/sqrt(1 - v²/c²).
pub fn lorentz_factor(v:f64)->f64{1.0/powers::sqrt(1.0-v*v/(units::SPEED_OF_LIGHT*units::SPEED_OF_LIGHT))}
/// Relativistic mass.
pub fn relativistic_mass(m0:f64,v:f64)->f64{m0*lorentz_factor(v)}
/// Mass-energy equivalence: E = mc².
pub fn mass_energy(m:f64)->f64{m*units::SPEED_OF_LIGHT*units::SPEED_OF_LIGHT}
/// Time dilation: t' = t*γ.
pub fn time_dilation(t:f64,v:f64)->f64{t*lorentz_factor(v)}
/// Length contraction: L' = L/γ.
pub fn length_contraction(l:f64,v:f64)->f64{l/lorentz_factor(v)}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_projectile(){assert!(approx_eq(projectile_range(10.0,PI/4.0,9.81),10.19,0.1));}
    #[test]fn test_kinetic(){assert!(approx_eq(kinetic_energy(2.0,3.0),9.0,1e-10));}
    #[test]fn test_ohm(){assert!(approx_eq(ohms_law_v(2.0,5.0),10.0,1e-10));}
    #[test]fn test_carnot(){assert!(approx_eq(carnot_efficiency(500.0,300.0),0.4,1e-10));}
    #[test]fn test_lorentz(){assert!(approx_eq(lorentz_factor(0.0),1.0,1e-10));}
    #[test]fn test_mass_energy(){let e=mass_energy(1.0);assert!(e>8e16&&e<9.1e16);}
}
