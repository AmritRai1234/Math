//! Bayesian inference: conjugate priors, MAP estimation, credible intervals.

use crate::powers;use crate::exponential;use crate::constants::PI;use crate::special;

/// Beta-Binomial conjugate update: posterior Beta(α', β') after observing k successes in n trials.
pub fn beta_binomial_update(prior_alpha:f64,prior_beta:f64,successes:u64,trials:u64)->(f64,f64){
    (prior_alpha+successes as f64,prior_beta+(trials-successes)as f64)
}

/// Posterior mean of Beta distribution.
pub fn beta_mean(alpha:f64,beta:f64)->f64{alpha/(alpha+beta)}

/// Posterior mode (MAP) of Beta distribution.
pub fn beta_mode(alpha:f64,beta:f64)->f64{if alpha>1.0&&beta>1.0{(alpha-1.0)/(alpha+beta-2.0)}else{beta_mean(alpha,beta)}}

/// Posterior variance of Beta distribution.
pub fn beta_variance(alpha:f64,beta:f64)->f64{alpha*beta/((alpha+beta)*(alpha+beta)*(alpha+beta+1.0))}

/// Normal-Normal conjugate: prior N(μ₀, σ₀²), data mean x̄ with n observations, known σ².
pub fn normal_normal_update(mu0:f64,sigma0_sq:f64,x_bar:f64,sigma_sq:f64,n:u64)->(f64,f64){
    let precision0=1.0/sigma0_sq;let precision_data=n as f64/sigma_sq;
    let new_precision=precision0+precision_data;
    let new_mean=(precision0*mu0+precision_data*x_bar)/new_precision;
    (new_mean,1.0/new_precision)
}

/// Gamma-Poisson conjugate: prior Gamma(α, β), observe n data points with sum s.
pub fn gamma_poisson_update(alpha:f64,beta:f64,n:u64,sum:f64)->(f64,f64){
    (alpha+sum,beta+n as f64)
}

/// Dirichlet-Multinomial conjugate update.
pub fn dirichlet_multinomial_update(prior_alpha:&[f64],counts:&[u64])->Vec<f64>{
    prior_alpha.iter().zip(counts.iter()).map(|(&a,&c)|a+c as f64).collect()
}

/// Dirichlet posterior mean.
pub fn dirichlet_mean(alpha:&[f64])->Vec<f64>{
    let total:f64=alpha.iter().sum();alpha.iter().map(|&a|a/total).collect()
}

/// Bayesian Information Criterion: BIC = k*ln(n) - 2*ln(L̂).
pub fn bic(log_likelihood:f64,k:u64,n:u64)->f64{
    k as f64*exponential::ln(n as f64)-2.0*log_likelihood
}

/// Akaike Information Criterion: AIC = 2k - 2*ln(L̂).
pub fn aic(log_likelihood:f64,k:u64)->f64{2.0*k as f64-2.0*log_likelihood}

/// Bayes factor: B₁₂ = P(data|M₁) / P(data|M₂).
pub fn bayes_factor(log_marginal_1:f64,log_marginal_2:f64)->f64{
    exponential::exp(log_marginal_1-log_marginal_2)
}

/// Jeffrey's prior for binomial: Beta(0.5, 0.5).
pub fn jeffreys_binomial()->(f64,f64){(0.5,0.5)}

/// Jeffrey's prior for Poisson: Gamma(0.5, 0).
pub fn jeffreys_poisson()->(f64,f64){(0.5,0.0)}

/// Highest Density Interval for Beta distribution (grid approximation).
pub fn beta_hdi(alpha:f64,beta_param:f64,credibility:f64)->(f64,f64){
    let n=1000;let dx=1.0/n as f64;
    let pdf=|x:f64|crate::probability::beta_pdf(x,alpha,beta_param);
    // Find the shortest interval containing credibility mass
    let mut best_width=1.0;let mut best_lo=0.0;let mut best_hi=1.0;
    for i in 0..n{let lo=i as f64*dx;let mut mass=0.0;
        for j in i..n{let hi=(j+1)as f64*dx;mass+=pdf((lo+hi)/2.0)*dx;
            if mass>=credibility{let width=hi-lo;if width<best_width{best_width=width;best_lo=lo;best_hi=hi;}break;}}}
    (best_lo,best_hi)
}

/// Posterior predictive for Beta-Binomial: P(X=k|data).
pub fn beta_binomial_predictive(k:u64,n:u64,alpha:f64,beta_param:f64)->f64{
    let binom=crate::combinatorics::combinations(n,k)as f64;
    let num=exponential::exp(special::lgamma(alpha+k as f64)+special::lgamma(beta_param+(n-k)as f64)-special::lgamma(alpha+beta_param+n as f64));
    let den=exponential::exp(special::lgamma(alpha)+special::lgamma(beta_param)-special::lgamma(alpha+beta_param));
    binom*num/den
}

/// Maximum A Posteriori estimate (grid search).
pub fn map_estimate(log_posterior:&dyn Fn(f64)->f64,lo:f64,hi:f64,n_grid:usize)->f64{
    let dx=(hi-lo)/n_grid as f64;let mut best_x=lo;let mut best_lp=f64::NEG_INFINITY;
    for i in 0..=n_grid{let x=lo+i as f64*dx;let lp=log_posterior(x);if lp>best_lp{best_lp=lp;best_x=x;}}best_x
}

/// Naive Bayes classifier: returns predicted class.
pub fn naive_bayes_predict(priors:&[f64],likelihoods:&[Vec<f64>])->usize{
    // likelihoods[class][feature] = P(feature|class)
    let posteriors:Vec<f64>=priors.iter().enumerate().map(|(c,&p)|{
        let ll:f64=likelihoods[c].iter().map(|&l|exponential::ln(l.max(1e-300))).sum();
        exponential::ln(p)+ll}).collect();
    posteriors.iter().enumerate().max_by(|(_,a),(_,b)|a.partial_cmp(b).unwrap()).map(|(i,_)|i).unwrap_or(0)
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_beta_update(){let(a,b)=beta_binomial_update(1.0,1.0,7,10);assert!(approx_eq(a,8.0,1e-10));assert!(approx_eq(b,4.0,1e-10));}
    #[test]fn test_beta_mean(){assert!(approx_eq(beta_mean(8.0,4.0),2.0/3.0,1e-6));}
    #[test]fn test_normal_update(){let(mu,_)=normal_normal_update(0.0,1.0,5.0,1.0,10);assert!(mu>4.0&&mu<5.5);}
    #[test]fn test_bic(){let b=bic(-100.0,3,100);assert!(b>0.0);}
}
