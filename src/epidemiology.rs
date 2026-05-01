//! Epidemiological models: SIR, SEIR, SIS.

/// SIR model: dS/dt = -βSI, dI/dt = βSI - γI, dR/dt = γI.
/// Returns Vec of (t, S, I, R) tuples.
pub fn sir(s0:f64,i0:f64,r0:f64,beta:f64,gamma:f64,dt:f64,steps:usize)->Vec<(f64,f64,f64,f64)>{
    let(mut s,mut i,mut r)=(s0,i0,r0);let n=s+i+r;let mut result=vec![(0.0,s,i,r)];
    for step in 1..=steps{let ds=-beta*s*i/n;let di=beta*s*i/n-gamma*i;let dr=gamma*i;
        s+=ds*dt;i+=di*dt;r+=dr*dt;result.push((step as f64*dt,s,i,r));}result
}

/// SEIR model with exposed compartment. Returns (t, S, E, I, R).
pub fn seir(s0:f64,e0:f64,i0:f64,r0:f64,beta:f64,sigma:f64,gamma:f64,dt:f64,steps:usize)->Vec<(f64,f64,f64,f64,f64)>{
    let(mut s,mut e,mut i,mut r)=(s0,e0,i0,r0);let n=s+e+i+r;let mut result=vec![(0.0,s,e,i,r)];
    for step in 1..=steps{let ds=-beta*s*i/n;let de=beta*s*i/n-sigma*e;let di=sigma*e-gamma*i;let dr=gamma*i;
        s+=ds*dt;e+=de*dt;i+=di*dt;r+=dr*dt;result.push((step as f64*dt,s,e,i,r));}result
}

/// SIS model (no immunity). Returns (t, S, I).
pub fn sis(s0:f64,i0:f64,beta:f64,gamma:f64,dt:f64,steps:usize)->Vec<(f64,f64,f64)>{
    let(mut s,mut i)=(s0,i0);let n=s+i;let mut result=vec![(0.0,s,i)];
    for step in 1..=steps{let ds=-beta*s*i/n+gamma*i;let di=beta*s*i/n-gamma*i;
        s+=ds*dt;i+=di*dt;result.push((step as f64*dt,s,i));}result
}

/// SIRD model (with deaths). Returns (t, S, I, R, D).
pub fn sird(s0:f64,i0:f64,r0:f64,d0:f64,beta:f64,gamma:f64,mu:f64,dt:f64,steps:usize)->Vec<(f64,f64,f64,f64,f64)>{
    let(mut s,mut i,mut r,mut d)=(s0,i0,r0,d0);let n=s+i+r;let mut result=vec![(0.0,s,i,r,d)];
    for step in 1..=steps{let ds=-beta*s*i/n;let di=beta*s*i/n-gamma*i-mu*i;let dr=gamma*i;let dd=mu*i;
        s+=ds*dt;i+=di*dt;r+=dr*dt;d+=dd*dt;result.push((step as f64*dt,s,i,r,d));}result
}

/// Basic reproduction number R0 = β/γ.
pub fn basic_reproduction_number(beta:f64,gamma:f64)->f64{beta/gamma}

/// Herd immunity threshold.
pub fn herd_immunity_threshold(r0:f64)->f64{1.0-1.0/r0}

/// Final epidemic size (transcendental equation: R∞ = 1 - exp(-R0 * R∞)).
pub fn final_size(r0:f64)->f64{
    let mut r=0.99;for _ in 0..1000{r=1.0-crate::exponential::exp(-r0*r);} r
}

/// Doubling time from growth rate.
pub fn doubling_time(growth_rate:f64)->f64{crate::exponential::ln(2.0)/growth_rate}

/// Effective reproduction number Re = R0 * S/N.
pub fn effective_r(r0:f64,susceptible:f64,total:f64)->f64{r0*susceptible/total}

/// Attack rate (fraction infected over epidemic).
pub fn attack_rate(total_infected:f64,population:f64)->f64{total_infected/population}

/// Case fatality rate.
pub fn case_fatality_rate(deaths:f64,cases:f64)->f64{deaths/cases}

/// Infection fatality rate.
pub fn infection_fatality_rate(deaths:f64,total_infected:f64)->f64{deaths/total_infected}

/// Serial interval exponential growth: I(t) = I0 * exp(r*t).
pub fn exponential_growth(i0:f64,rate:f64,t:f64)->f64{i0*crate::exponential::exp(rate*t)}

/// Logistic growth: N(t) = K / (1 + ((K-N0)/N0) * exp(-r*t)).
pub fn logistic_growth(n0:f64,k:f64,r:f64,t:f64)->f64{k/(1.0+((k-n0)/n0)*crate::exponential::exp(-r*t))}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_sir(){let result=sir(999.0,1.0,0.0,0.3,0.1,0.1,1000);
        let last=result.last().unwrap();assert!(approx_eq(last.1+last.2+last.3,1000.0,1.0));}
    #[test]fn test_r0(){assert!(approx_eq(basic_reproduction_number(0.3,0.1),3.0,1e-10));}
    #[test]fn test_herd(){assert!(approx_eq(herd_immunity_threshold(3.0),2.0/3.0,1e-10));}
    #[test]fn test_doubling(){assert!(approx_eq(doubling_time(0.1),crate::exponential::ln(2.0)/0.1,1e-10));}
}
