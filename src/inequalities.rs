//! Mathematical inequalities: AM-GM, Cauchy-Schwarz, Jensen, power means, triangle.

use crate::powers;use crate::arithmetic::abs;

/// Arithmetic mean.
pub fn arithmetic_mean(data:&[f64])->f64{data.iter().sum::<f64>()/data.len()as f64}

/// Geometric mean.
pub fn geometric_mean(data:&[f64])->f64{
    let n=data.len()as f64;let log_sum:f64=data.iter().map(|&x|crate::exponential::ln(x)).sum();
    crate::exponential::exp(log_sum/n)
}

/// Harmonic mean.
pub fn harmonic_mean(data:&[f64])->f64{
    let n=data.len()as f64;n/data.iter().map(|&x|1.0/x).sum::<f64>()
}

/// Quadratic mean (RMS).
pub fn quadratic_mean(data:&[f64])->f64{powers::sqrt(data.iter().map(|&x|x*x).sum::<f64>()/data.len()as f64)}

/// Power mean M_p.
pub fn power_mean(data:&[f64],p:f64)->f64{
    if p==0.0{return geometric_mean(data);}
    powers::pow(data.iter().map(|&x|powers::pow(x,p)).sum::<f64>()/data.len()as f64,1.0/p)
}

/// Verify AM-GM inequality: AM ≥ GM (for positive numbers).
pub fn verify_am_gm(data:&[f64])->bool{arithmetic_mean(data)>=geometric_mean(data)-1e-10}

/// Verify HM ≤ GM ≤ AM ≤ QM chain.
pub fn verify_mean_inequality(data:&[f64])->bool{
    let hm=harmonic_mean(data);let gm=geometric_mean(data);
    let am=arithmetic_mean(data);let qm=quadratic_mean(data);
    hm<=gm+1e-10&&gm<=am+1e-10&&am<=qm+1e-10
}

/// Cauchy-Schwarz inequality: (Σa_i*b_i)² ≤ (Σa_i²)(Σb_i²).
pub fn cauchy_schwarz_lhs(a:&[f64],b:&[f64])->f64{
    let dot:f64=a.iter().zip(b.iter()).map(|(&x,&y)|x*y).sum();dot*dot
}
pub fn cauchy_schwarz_rhs(a:&[f64],b:&[f64])->f64{
    let sa:f64=a.iter().map(|&x|x*x).sum();let sb:f64=b.iter().map(|&x|x*x).sum();sa*sb
}
pub fn verify_cauchy_schwarz(a:&[f64],b:&[f64])->bool{cauchy_schwarz_lhs(a,b)<=cauchy_schwarz_rhs(a,b)+1e-10}

/// Triangle inequality: |a + b| ≤ |a| + |b| (for vectors).
pub fn verify_triangle_inequality(a:&[f64],b:&[f64])->bool{
    let sum_norm:f64=powers::sqrt(a.iter().zip(b.iter()).map(|(&x,&y)|(x+y)*(x+y)).sum());
    let norm_a:f64=powers::sqrt(a.iter().map(|&x|x*x).sum());
    let norm_b:f64=powers::sqrt(b.iter().map(|&x|x*x).sum());
    sum_norm<=norm_a+norm_b+1e-10
}

/// Minkowski inequality: ||a+b||_p ≤ ||a||_p + ||b||_p.
pub fn lp_norm(data:&[f64],p:f64)->f64{
    powers::pow(data.iter().map(|&x|powers::pow(abs(x),p)).sum::<f64>(),1.0/p)
}
pub fn verify_minkowski(a:&[f64],b:&[f64],p:f64)->bool{
    let sum:Vec<f64>=a.iter().zip(b.iter()).map(|(&x,&y)|x+y).collect();
    lp_norm(&sum,p)<=lp_norm(a,p)+lp_norm(b,p)+1e-10
}

/// Hölder's inequality: Σ|a_i*b_i| ≤ ||a||_p * ||b||_q where 1/p + 1/q = 1.
pub fn verify_holder(a:&[f64],b:&[f64],p:f64)->bool{
    let q=p/(p-1.0);
    let lhs:f64=a.iter().zip(b.iter()).map(|(&x,&y)|abs(x*y)).sum();
    lhs<=lp_norm(a,p)*lp_norm(b,q)+1e-10
}

/// Jensen's inequality for convex function: f(E[X]) ≤ E[f(X)].
pub fn verify_jensen_convex(f:&dyn Fn(f64)->f64,data:&[f64])->bool{
    let mean=arithmetic_mean(data);let f_mean=f(mean);
    let mean_f=arithmetic_mean(&data.iter().map(|&x|f(x)).collect::<Vec<_>>());
    f_mean<=mean_f+1e-10
}

/// Chebyshev's sum inequality: if a and b are similarly sorted,
/// n*Σa_i*b_i ≥ (Σa_i)*(Σb_i).
pub fn verify_chebyshev_sum(a:&[f64],b:&[f64])->bool{
    let n=a.len()as f64;let dot:f64=a.iter().zip(b.iter()).map(|(&x,&y)|x*y).sum();
    let sa:f64=a.iter().sum();let sb:f64=b.iter().sum();
    n*dot>=sa*sb-1e-10
}

/// Weighted power mean.
pub fn weighted_power_mean(data:&[f64],weights:&[f64],p:f64)->f64{
    let total_w:f64=weights.iter().sum();
    if p==0.0{// Weighted geometric mean
        let log_sum:f64=data.iter().zip(weights.iter()).map(|(&x,&w)|w*crate::exponential::ln(x)).sum();
        return crate::exponential::exp(log_sum/total_w);
    }
    powers::pow(data.iter().zip(weights.iter()).map(|(&x,&w)|w*powers::pow(x,p)).sum::<f64>()/total_w,1.0/p)
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_am_gm(){assert!(verify_am_gm(&[2.0,8.0]));}
    #[test]fn test_cauchy_schwarz(){assert!(verify_cauchy_schwarz(&[1.0,2.0,3.0],&[4.0,5.0,6.0]));}
    #[test]fn test_triangle(){assert!(verify_triangle_inequality(&[1.0,2.0],&[3.0,4.0]));}
    #[test]fn test_mean_chain(){assert!(verify_mean_inequality(&[1.0,2.0,3.0,4.0]));}
    #[test]fn test_power_mean(){assert!(approx_eq(power_mean(&[4.0,4.0],2.0),4.0,1e-10));}
    #[test]fn test_jensen(){assert!(verify_jensen_convex(&|x|x*x,&[1.0,2.0,3.0,4.0]));}
}
