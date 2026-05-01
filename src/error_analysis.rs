//! Numerical error analysis and propagation.

use crate::arithmetic::abs;
use crate::powers;

/// Machine epsilon for f64.
pub const MACHINE_EPS:f64=f64::EPSILON; // ≈ 2.22e-16

/// Absolute error.
pub fn absolute_error(exact:f64,approx:f64)->f64{abs(exact-approx)}

/// Relative error.
pub fn relative_error(exact:f64,approx:f64)->f64{if exact==0.0{abs(approx)}else{abs((exact-approx)/exact)}}

/// Ulps (units in the last place) between two f64 values.
pub fn ulps_distance(a:f64,b:f64)->u64{
    let ai=a.to_bits()as i64;let bi=b.to_bits()as i64;(ai-bi).unsigned_abs()
}

/// Error propagation for f(x) ≈ f(x0) + f'(x0)*δx.
/// Returns propagated error given input uncertainties.
pub fn propagate_error(partials:&[f64],uncertainties:&[f64])->f64{
    powers::sqrt(partials.iter().zip(uncertainties.iter()).map(|(&p,&u)|p*p*u*u).sum())
}

/// Error propagation for addition: δ(a+b) = sqrt(δa² + δb²).
pub fn error_add(da:f64,db:f64)->f64{powers::sqrt(da*da+db*db)}

/// Error propagation for multiplication: δ(a*b)/|a*b| = sqrt((δa/a)² + (δb/b)²).
pub fn error_mul(a:f64,da:f64,b:f64,db:f64)->f64{
    abs(a*b)*powers::sqrt((da/a)*(da/a)+(db/b)*(db/b))
}

/// Error propagation for division.
pub fn error_div(a:f64,da:f64,b:f64,db:f64)->f64{
    abs(a/b)*powers::sqrt((da/a)*(da/a)+(db/b)*(db/b))
}

/// Error propagation for power: δ(a^n) = |n * a^(n-1)| * δa.
pub fn error_pow(a:f64,da:f64,n:f64)->f64{abs(n*powers::pow(a,n-1.0))*da}

/// Condition number of a function at point x: |x*f'(x)/f(x)|.
pub fn condition_number(f:&dyn Fn(f64)->f64,x:f64)->f64{
    let h=1e-8*abs(x).max(1e-8);let fp=(f(x+h)-f(x-h))/(2.0*h);
    let fx=f(x);if fx==0.0{return f64::INFINITY;}abs(x*fp/fx)
}

/// Numerical stability analysis: compute forward and backward error.
pub fn stability_analysis(exact:f64,computed:f64,input:f64,perturbed_input:f64)->(f64,f64){
    let forward=relative_error(exact,computed);
    let backward=relative_error(input,perturbed_input);
    (forward,backward)
}

/// Kahan summation (compensated summation for better accuracy).
pub fn kahan_sum(data:&[f64])->f64{
    let mut sum=0.0;let mut c=0.0; // compensation
    for &x in data{let y=x-c;let t=sum+y;c=(t-sum)-y;sum=t;}sum
}

/// Pairwise summation (better than naive, simpler than Kahan).
pub fn pairwise_sum(data:&[f64])->f64{
    if data.len()<=128{return data.iter().sum();}
    let mid=data.len()/2;pairwise_sum(&data[..mid])+pairwise_sum(&data[mid..])
}

/// Significant digits of agreement between two values.
pub fn significant_digits(a:f64,b:f64)->u32{
    if a==b{return 16;}let re=relative_error(a,b);if re==0.0{return 16;}
    let digits=-crate::exponential::log10(re);if digits<0.0{0}else{digits as u32}
}

/// Rounding error bound for n additions.
pub fn summation_error_bound(n:usize)->f64{n as f64*MACHINE_EPS/(1.0-n as f64*MACHINE_EPS)}

/// Check if two f64 values are equal within n ULPs.
pub fn ulps_eq(a:f64,b:f64,max_ulps:u64)->bool{ulps_distance(a,b)<=max_ulps}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_absolute(){assert!(approx_eq(absolute_error(1.0,1.001),0.001,1e-10));}
    #[test]fn test_relative(){assert!(approx_eq(relative_error(100.0,99.0),0.01,1e-10));}
    #[test]fn test_propagate(){let p=[1.0,1.0];let u=[0.1,0.2];
        assert!(approx_eq(propagate_error(&p,&u),powers::sqrt(0.05),1e-10));}
    #[test]fn test_kahan(){let data:Vec<f64>=(0..10000).map(|_|0.1).collect();
        let naive:f64=data.iter().sum();let kahan=kahan_sum(&data);
        assert!(abs(kahan-1000.0)<abs(naive-1000.0)||abs(kahan-1000.0)<1e-10);}
    #[test]fn test_ulps(){assert!(ulps_eq(1.0,1.0+f64::EPSILON,1));}
}
