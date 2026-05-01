//! Stochastic calculus: Brownian motion, SDEs, Monte Carlo.

use crate::powers;
use crate::exponential;

/// Generate Brownian motion path (Wiener process).
pub fn brownian_motion(n:usize,dt:f64,seed:u64)->Vec<f64>{
    let mut rng=crate::random::Rng::new(seed);let sqrt_dt=powers::sqrt(dt);
    let mut path=vec![0.0];let mut val=0.0;
    for _ in 0..n{let(z,_)=crate::random::box_muller(&mut rng);val+=sqrt_dt*z;path.push(val);}path
}

/// Geometric Brownian motion: S(t) = S0 * exp((μ - σ²/2)t + σW(t)).
pub fn geometric_brownian(s0:f64,mu:f64,sigma:f64,dt:f64,n:usize,seed:u64)->Vec<f64>{
    let mut rng=crate::random::Rng::new(seed);let sqrt_dt=powers::sqrt(dt);
    let mut path=vec![s0];let mut s=s0;
    for _ in 0..n{let(z,_)=crate::random::box_muller(&mut rng);
        s*=exponential::exp((mu-0.5*sigma*sigma)*dt+sigma*sqrt_dt*z);path.push(s);}path
}

/// Ornstein-Uhlenbeck process: dX = θ(μ - X)dt + σdW.
pub fn ornstein_uhlenbeck(theta:f64,mu:f64,sigma:f64,x0:f64,dt:f64,n:usize,seed:u64)->Vec<f64>{
    let mut rng=crate::random::Rng::new(seed);let sqrt_dt=powers::sqrt(dt);
    let mut path=vec![x0];let mut x=x0;
    for _ in 0..n{let(z,_)=crate::random::box_muller(&mut rng);
        x+=theta*(mu-x)*dt+sigma*sqrt_dt*z;path.push(x);}path
}

/// Euler-Maruyama SDE solver: dX = a(X,t)dt + b(X,t)dW.
pub fn euler_maruyama(a:&dyn Fn(f64,f64)->f64,b:&dyn Fn(f64,f64)->f64,x0:f64,t0:f64,dt:f64,n:usize,seed:u64)->Vec<(f64,f64)>{
    let mut rng=crate::random::Rng::new(seed);let sqrt_dt=powers::sqrt(dt);
    let mut path=vec![(t0,x0)];let mut x=x0;let mut t=t0;
    for _ in 0..n{let(z,_)=crate::random::box_muller(&mut rng);
        x+=a(x,t)*dt+b(x,t)*sqrt_dt*z;t+=dt;path.push((t,x));}path
}

/// Milstein method (higher order SDE solver).
pub fn milstein(a:&dyn Fn(f64,f64)->f64,b:&dyn Fn(f64,f64)->f64,db:&dyn Fn(f64,f64)->f64,x0:f64,t0:f64,dt:f64,n:usize,seed:u64)->Vec<(f64,f64)>{
    let mut rng=crate::random::Rng::new(seed);let sqrt_dt=powers::sqrt(dt);
    let mut path=vec![(t0,x0)];let mut x=x0;let mut t=t0;
    for _ in 0..n{let(z,_)=crate::random::box_muller(&mut rng);let dw=sqrt_dt*z;
        x+=a(x,t)*dt+b(x,t)*dw+0.5*b(x,t)*db(x,t)*(dw*dw-dt);t+=dt;path.push((t,x));}path
}

/// Monte Carlo option pricing (European call).
pub fn monte_carlo_call(s0:f64,k:f64,r:f64,sigma:f64,t:f64,n_paths:usize,n_steps:usize,seed:u64)->f64{
    let dt=t/n_steps as f64;let mut sum=0.0;
    for i in 0..n_paths{let path=geometric_brownian(s0,r,sigma,dt,n_steps,seed.wrapping_add(i as u64));
        let st=*path.last().unwrap();let payoff=if st>k{st-k}else{0.0};sum+=payoff;}
    exponential::exp(-r*t)*sum/n_paths as f64
}

/// Monte Carlo option pricing (European put).
pub fn monte_carlo_put(s0:f64,k:f64,r:f64,sigma:f64,t:f64,n_paths:usize,n_steps:usize,seed:u64)->f64{
    let dt=t/n_steps as f64;let mut sum=0.0;
    for i in 0..n_paths{let path=geometric_brownian(s0,r,sigma,dt,n_steps,seed.wrapping_add(i as u64));
        let st=*path.last().unwrap();let payoff=if k>st{k-st}else{0.0};sum+=payoff;}
    exponential::exp(-r*t)*sum/n_paths as f64
}

/// Black-Scholes formula (European call).
pub fn black_scholes_call(s:f64,k:f64,r:f64,sigma:f64,t:f64)->f64{
    let d1=(exponential::ln(s/k)+(r+0.5*sigma*sigma)*t)/(sigma*powers::sqrt(t));
    let d2=d1-sigma*powers::sqrt(t);
    s*normal_cdf(d1)-k*exponential::exp(-r*t)*normal_cdf(d2)
}

/// Black-Scholes formula (European put).
pub fn black_scholes_put(s:f64,k:f64,r:f64,sigma:f64,t:f64)->f64{
    let d1=(exponential::ln(s/k)+(r+0.5*sigma*sigma)*t)/(sigma*powers::sqrt(t));
    let d2=d1-sigma*powers::sqrt(t);
    k*exponential::exp(-r*t)*normal_cdf(-d2)-s*normal_cdf(-d1)
}

fn normal_cdf(x:f64)->f64{0.5*(1.0+crate::special::erf(x/crate::constants::SQRT2))}

/// Variance of Monte Carlo estimator.
pub fn mc_variance(samples:&[f64])->f64{
    let n=samples.len()as f64;let mean=samples.iter().sum::<f64>()/n;
    samples.iter().map(|&x|(x-mean)*(x-mean)).sum::<f64>()/(n-1.0)/n
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_brownian(){let path=brownian_motion(1000,0.01,42);assert_eq!(path.len(),1001);assert_eq!(path[0],0.0);}
    #[test]fn test_gbm(){let path=geometric_brownian(100.0,0.05,0.2,0.01,100,42);assert!(path[0]==100.0);assert!(path.last().unwrap()>&0.0);}
    #[test]fn test_black_scholes(){let c=black_scholes_call(100.0,100.0,0.05,0.2,1.0);assert!(c>5.0&&c<15.0);}
    #[test]fn test_put_call_parity(){let s=100.0;let k=100.0;let r=0.05;let sigma=0.2;let t=1.0;
        let c=black_scholes_call(s,k,r,sigma,t);let p=black_scholes_put(s,k,r,sigma,t);
        let parity=c-p-(s-k*exponential::exp(-r*t));assert!(parity.abs()<1e-6);}
}
