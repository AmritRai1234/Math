//! Queueing theory: M/M/1, M/M/c, Erlang formulas, Little's law.

use crate::combinatorics::factorial;
use crate::powers;

/// M/M/1 queue metrics. λ = arrival rate, μ = service rate.
pub fn mm1(lambda:f64,mu:f64)->Mm1Result{
    let rho=lambda/mu;
    Mm1Result{rho,l_q:rho*rho/(1.0-rho),l_s:rho/(1.0-rho),
        w_q:rho/(mu-lambda),w_s:1.0/(mu-lambda),p0:1.0-rho}
}

#[derive(Debug,Clone)]
pub struct Mm1Result{pub rho:f64,pub l_q:f64,pub l_s:f64,pub w_q:f64,pub w_s:f64,pub p0:f64}

/// M/M/c queue metrics.
pub fn mmc(lambda:f64,mu:f64,c:u64)->MmcResult{
    let rho=lambda/(c as f64*mu);let a=lambda/mu;
    let sum:f64=(0..c).map(|k|powers::pow(a,k as f64)/factorial(k)as f64).sum();
    let last_term=powers::pow(a,c as f64)/(factorial(c)as f64*(1.0-rho));
    let p0=1.0/(sum+last_term);
    let c_prob=powers::pow(a,c as f64)*p0/(factorial(c)as f64*(1.0-rho));
    let l_q=c_prob*rho/(1.0-rho);let w_q=l_q/lambda;let w_s=w_q+1.0/mu;let l_s=lambda*w_s;
    MmcResult{rho,l_q,l_s,w_q,w_s,p0,c_prob}
}

#[derive(Debug,Clone)]
pub struct MmcResult{pub rho:f64,pub l_q:f64,pub l_s:f64,pub w_q:f64,pub w_s:f64,pub p0:f64,pub c_prob:f64}

/// Little's law: L = λ * W.
pub fn little_l(lambda:f64,w:f64)->f64{lambda*w}
pub fn little_w(lambda:f64,l:f64)->f64{l/lambda}

/// Erlang B formula: probability of blocking with c servers, traffic intensity a.
pub fn erlang_b(a:f64,c:u64)->f64{
    let num=powers::pow(a,c as f64)/factorial(c)as f64;
    let den:f64=(0..=c).map(|k|powers::pow(a,k as f64)/factorial(k)as f64).sum();
    num/den
}

/// Erlang C formula: probability of waiting with c servers.
pub fn erlang_c(a:f64,c:u64)->f64{
    let rho=a/c as f64;if rho>=1.0{return 1.0;}
    let num=powers::pow(a,c as f64)/(factorial(c)as f64*(1.0-rho));
    let sum:f64=(0..c).map(|k|powers::pow(a,k as f64)/factorial(k)as f64).sum();
    num/(sum+num)
}

/// M/M/1/K (finite buffer). P(n customers) for n=0..K.
pub fn mm1k_probabilities(lambda:f64,mu:f64,k:u64)->Vec<f64>{
    let rho=lambda/mu;
    if(rho-1.0).abs()<1e-10{return vec![1.0/(k+1)as f64;(k+1)as usize];}
    let p0=(1.0-rho)/(1.0-powers::pow(rho,(k+1)as f64));
    (0..=k).map(|n|p0*powers::pow(rho,n as f64)).collect()
}

/// M/M/∞ queue (infinite servers). P(n) = Poisson with mean λ/μ.
pub fn mm_inf_prob(lambda:f64,mu:f64,n:u64)->f64{
    let a=lambda/mu;crate::distributions::poisson_pmf(n,a)
}

/// Utilization factor.
pub fn utilization(lambda:f64,mu:f64,c:u64)->f64{lambda/(c as f64*mu)}

/// Average service time.
pub fn avg_service_time(mu:f64)->f64{1.0/mu}

/// Throughput (effective arrival rate for finite systems).
pub fn throughput(lambda:f64,p_loss:f64)->f64{lambda*(1.0-p_loss)}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_mm1(){let r=mm1(2.0,3.0);assert!(approx_eq(r.rho,2.0/3.0,1e-10));assert!(approx_eq(r.l_s,2.0,1e-10));}
    #[test]fn test_little(){assert!(approx_eq(little_l(5.0,2.0),10.0,1e-10));}
    #[test]fn test_erlang_b(){let pb=erlang_b(10.0,12);assert!(pb>0.0&&pb<1.0);}
    #[test]fn test_utilization(){assert!(approx_eq(utilization(4.0,3.0,2),2.0/3.0,1e-10));}
}
