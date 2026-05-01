//! Probability theory: Bayes, conditional, joint, total probability, copulas, stochastic dominance.

use crate::powers;use crate::exponential;use crate::combinatorics;use crate::constants::PI;

// ─── Fundamental Rules ───

/// Bayes' theorem: P(A|B) = P(B|A) * P(A) / P(B).
pub fn bayes(p_b_given_a:f64,p_a:f64,p_b:f64)->f64{p_b_given_a*p_a/p_b}

/// Law of total probability: P(B) = Σ P(B|Ai) * P(Ai).
pub fn total_probability(p_b_given_a:&[f64],p_a:&[f64])->f64{
    p_b_given_a.iter().zip(p_a.iter()).map(|(&b,&a)|b*a).sum()
}

/// Chain rule: P(A∩B∩C) = P(A) * P(B|A) * P(C|A∩B).
pub fn chain_rule(probs:&[f64])->f64{probs.iter().product()}

/// Inclusion-exclusion for 2 events: P(A∪B) = P(A) + P(B) - P(A∩B).
pub fn inclusion_exclusion_2(pa:f64,pb:f64,p_ab:f64)->f64{pa+pb-p_ab}

/// Inclusion-exclusion for 3 events.
pub fn inclusion_exclusion_3(pa:f64,pb:f64,pc:f64,pab:f64,pac:f64,pbc:f64,pabc:f64)->f64{
    pa+pb+pc-pab-pac-pbc+pabc
}

/// P(A|B) from joint and marginal.
pub fn conditional(p_ab:f64,p_b:f64)->f64{if p_b>0.0{p_ab/p_b}else{0.0}}

/// Independence check: P(A∩B) ≈ P(A)*P(B)?
pub fn are_independent(p_a:f64,p_b:f64,p_ab:f64,tol:f64)->bool{(p_ab-p_a*p_b).abs()<tol}

// ─── Combinatorial Probability ───

/// Birthday problem: probability that at least 2 share a birthday in n people.
pub fn birthday_problem(n:u32)->f64{
    let mut p_no_match=1.0;for i in 1..n{p_no_match*=(365.0-i as f64)/365.0;}1.0-p_no_match
}

/// Coupon collector: expected number of draws to collect all n types.
pub fn coupon_collector_expected(n:u64)->f64{
    let mut sum=0.0;for i in 1..=n{sum+=1.0/i as f64;}n as f64*sum
}

/// Derangement probability: P(no fixed points) = Σ(-1)^k/k!
pub fn derangement_probability(n:u64)->f64{
    (0..=n).map(|k|{let sign=if k%2==0{1.0}else{-1.0};sign/combinatorics::factorial(k)as f64}).sum()
}

/// Runs test: expected number of runs in sequence of n1 ones and n2 zeros.
pub fn expected_runs(n1:u64,n2:u64)->f64{1.0+2.0*n1 as f64*n2 as f64/(n1+n2)as f64}

// ─── Distributions & Density Functions ───

/// Normal PDF.
pub fn normal_pdf(x:f64,mu:f64,sigma:f64)->f64{
    exponential::exp(-0.5*((x-mu)/sigma)*((x-mu)/sigma))/(sigma*powers::sqrt(2.0*PI))
}

/// Normal CDF (using erf).
pub fn normal_cdf(x:f64,mu:f64,sigma:f64)->f64{
    0.5*(1.0+crate::special::erf((x-mu)/(sigma*crate::constants::SQRT2)))
}

/// Inverse normal CDF (probit function, rational approximation).
pub fn normal_quantile(p:f64)->f64{
    if p<=0.0{return f64::NEG_INFINITY;}if p>=1.0{return f64::INFINITY;}
    // Rational approximation (Abramowitz & Stegun)
    let t=if p<0.5{powers::sqrt(-2.0*exponential::ln(p))}else{powers::sqrt(-2.0*exponential::ln(1.0-p))};
    let c=[2.515517,0.802853,0.010328];let d=[1.432788,0.189269,0.001308];
    let result=t-(c[0]+c[1]*t+c[2]*t*t)/(1.0+d[0]*t+d[1]*t*t+d[2]*t*t*t);
    if p<0.5{-result}else{result}
}

/// Log-normal PDF.
pub fn lognormal_pdf(x:f64,mu:f64,sigma:f64)->f64{
    if x<=0.0{return 0.0;}
    exponential::exp(-0.5*((exponential::ln(x)-mu)/sigma)*((exponential::ln(x)-mu)/sigma))/(x*sigma*powers::sqrt(2.0*PI))
}

/// Exponential PDF.
pub fn exponential_pdf(x:f64,lambda:f64)->f64{if x<0.0{0.0}else{lambda*exponential::exp(-lambda*x)}}
/// Exponential CDF.
pub fn exponential_cdf(x:f64,lambda:f64)->f64{if x<0.0{0.0}else{1.0-exponential::exp(-lambda*x)}}

/// Gamma PDF (using Stirling's approximation for Γ).
pub fn gamma_pdf(x:f64,k:f64,theta:f64)->f64{
    if x<=0.0{return 0.0;}
    let log_gamma_k=crate::special::lgamma(k);
    exponential::exp((k-1.0)*exponential::ln(x)-x/theta-k*exponential::ln(theta)-log_gamma_k)
}

/// Beta PDF.
pub fn beta_pdf(x:f64,alpha:f64,beta:f64)->f64{
    if x<0.0||x>1.0{return 0.0;}
    let b=crate::special::lgamma(alpha)+crate::special::lgamma(beta)-crate::special::lgamma(alpha+beta);
    exponential::exp((alpha-1.0)*exponential::ln(x.max(1e-300))+(beta-1.0)*exponential::ln((1.0-x).max(1e-300))-b)
}

/// Weibull PDF.
pub fn weibull_pdf(x:f64,k:f64,lambda:f64)->f64{
    if x<0.0{return 0.0;}(k/lambda)*powers::pow(x/lambda,k-1.0)*exponential::exp(-powers::pow(x/lambda,k))
}

/// Cauchy PDF.
pub fn cauchy_pdf(x:f64,x0:f64,gamma:f64)->f64{
    1.0/(PI*gamma*(1.0+((x-x0)/gamma)*((x-x0)/gamma)))
}

/// Pareto PDF.
pub fn pareto_pdf(x:f64,x_m:f64,alpha:f64)->f64{
    if x<x_m{0.0}else{alpha*powers::pow(x_m,alpha)/powers::pow(x,alpha+1.0)}
}

/// Student's t PDF.
pub fn student_t_pdf(x:f64,nu:f64)->f64{
    let coeff=exponential::exp(crate::special::lgamma((nu+1.0)/2.0)-crate::special::lgamma(nu/2.0))
        /(powers::sqrt(nu*PI));
    coeff*powers::pow(1.0+x*x/nu,-(nu+1.0)/2.0)
}

/// Chi-squared PDF.
pub fn chi_squared_pdf(x:f64,k:f64)->f64{gamma_pdf(x,k/2.0,2.0)}

/// F-distribution PDF.
pub fn f_pdf(x:f64,d1:f64,d2:f64)->f64{
    if x<=0.0{return 0.0;}
    let num=powers::pow(d1*x,d1)*powers::pow(d2,d2)/powers::pow(d1*x+d2,d1+d2);
    let den=x*exponential::exp(crate::special::lgamma(d1/2.0)+crate::special::lgamma(d2/2.0)-crate::special::lgamma((d1+d2)/2.0));
    if den.abs()<1e-300{0.0}else{powers::sqrt(num/den.abs())}
}

// ─── Probability Inequalities ───

/// Markov's inequality: P(X ≥ a) ≤ E[X]/a for X ≥ 0.
pub fn markov_bound(expected_value:f64,a:f64)->f64{(expected_value/a).min(1.0)}

/// Chebyshev's inequality: P(|X-μ| ≥ kσ) ≤ 1/k².
pub fn chebyshev_bound(k:f64)->f64{(1.0/(k*k)).min(1.0)}

/// Chernoff bound (multiplicative, for sum of Bernoullis): P(X ≥ (1+δ)μ).
pub fn chernoff_bound(mu:f64,delta:f64)->f64{
    exponential::exp(-mu*delta*delta/3.0).min(1.0)
}

/// Hoeffding's inequality: P(|X̄ - E[X̄]| ≥ t) ≤ 2*exp(-2nt²/(b-a)²).
pub fn hoeffding_bound(n:usize,t:f64,a:f64,b:f64)->f64{
    (2.0*exponential::exp(-2.0*n as f64*t*t/((b-a)*(b-a)))).min(1.0)
}

// ─── Moments & Generating Functions ───

/// Raw moment E[X^k] from samples.
pub fn raw_moment(data:&[f64],k:u32)->f64{data.iter().map(|&x|powers::pow(x,k as f64)).sum::<f64>()/data.len()as f64}

/// Central moment E[(X-μ)^k] from samples.
pub fn central_moment(data:&[f64],k:u32)->f64{
    let mu=crate::stats::mean(data);data.iter().map(|&x|powers::pow(x-mu,k as f64)).sum::<f64>()/data.len()as f64
}

/// Moment generating function of normal at t: M(t) = exp(μt + σ²t²/2).
pub fn mgf_normal(t:f64,mu:f64,sigma:f64)->f64{exponential::exp(mu*t+0.5*sigma*sigma*t*t)}

/// Characteristic function of normal at t: φ(t) = exp(iμt - σ²t²/2).
pub fn cf_normal_mag(t:f64,sigma:f64)->f64{exponential::exp(-0.5*sigma*sigma*t*t)}

// ─── Copulas ───

/// Gaussian copula (bivariate): C(u,v) with correlation ρ.
pub fn gaussian_copula(u:f64,v:f64,rho:f64)->f64{
    let x=normal_quantile(u);let y=normal_quantile(v);
    normal_cdf((x-rho*y)/powers::sqrt(1.0-rho*rho),0.0,1.0)
}

/// Clayton copula: C(u,v) = (u^{-θ} + v^{-θ} - 1)^{-1/θ}.
pub fn clayton_copula(u:f64,v:f64,theta:f64)->f64{
    powers::pow(powers::pow(u,-theta)+powers::pow(v,-theta)-1.0,-1.0/theta)
}

// ─── Stochastic Orders ───

/// First order stochastic dominance check (empirical CDFs).
pub fn first_order_dominance(a:&[f64],b:&[f64])->bool{
    let mut sa=a.to_vec();sa.sort_by(|x,y|x.partial_cmp(y).unwrap());
    let mut sb=b.to_vec();sb.sort_by(|x,y|x.partial_cmp(y).unwrap());
    let n=sa.len().max(sb.len());
    (0..n).all(|i|{let pa=(i+1)as f64/sa.len()as f64;let pb=(i+1)as f64/sb.len()as f64;pa<=pb})
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_bayes(){// P(disease|positive) with P(pos|disease)=0.99, P(disease)=0.01, P(pos)=0.0297
        let p=bayes(0.99,0.01,0.0297);assert!(p>0.33&&p<0.34);}
    #[test]fn test_birthday(){assert!(birthday_problem(23)>0.50);}
    #[test]fn test_normal_pdf(){assert!(approx_eq(normal_pdf(0.0,0.0,1.0),0.3989422804,1e-6));}
    #[test]fn test_normal_cdf(){assert!(approx_eq(normal_cdf(0.0,0.0,1.0),0.5,1e-6));}
    #[test]fn test_chebyshev(){assert!(approx_eq(chebyshev_bound(2.0),0.25,1e-10));}
    #[test]fn test_derangement(){let p=derangement_probability(10);assert!((p-1.0/std::f64::consts::E).abs()<0.001);}
    #[test]fn test_total_prob(){let p=total_probability(&[0.8,0.2],&[0.5,0.5]);assert!(approx_eq(p,0.5,1e-10));}
}
