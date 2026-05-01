//! Regression and curve fitting.

use crate::powers;
use crate::stats;

/// Simple linear regression: y = a + b*x. Returns (intercept, slope).
pub fn linear_regression(x:&[f64],y:&[f64])->(f64,f64){
    let n=x.len()as f64;let mx=stats::mean(x);let my=stats::mean(y);
    let num:f64=x.iter().zip(y.iter()).map(|(&xi,&yi)|(xi-mx)*(yi-my)).sum();
    let den:f64=x.iter().map(|&xi|(xi-mx)*(xi-mx)).sum();
    let b=num/den;let a=my-b*mx;(a,b)
}

/// R² (coefficient of determination).
pub fn r_squared(x:&[f64],y:&[f64])->f64{
    let(a,b)=linear_regression(x,y);let my=stats::mean(y);
    let ss_res:f64=x.iter().zip(y.iter()).map(|(&xi,&yi)|{let pred=a+b*xi;(yi-pred)*(yi-pred)}).sum();
    let ss_tot:f64=y.iter().map(|&yi|(yi-my)*(yi-my)).sum();
    1.0-ss_res/ss_tot
}

/// Polynomial regression of degree d. Returns coefficients [a0, a1, ..., ad].
pub fn polynomial_regression(x:&[f64],y:&[f64],degree:usize)->Vec<f64>{
    let n=x.len();let m=degree+1;
    // Build Vandermonde matrix X^T*X and X^T*y
    let mut xtx=vec![vec![0.0;m];m];let mut xty=vec![0.0;m];
    for i in 0..m{for j in 0..m{for k in 0..n{xtx[i][j]+=powers::pow(x[k],(i+j)as f64);}}
        for k in 0..n{xty[i]+=powers::pow(x[k],i as f64)*y[k];}}
    // Solve via Gaussian elimination
    crate::linalg::solve_linear(&xtx,&xty)
}

/// Evaluate polynomial regression at x.
pub fn poly_predict(coeffs:&[f64],x:f64)->f64{
    coeffs.iter().enumerate().map(|(i,&c)|c*powers::pow(x,i as f64)).sum()
}

/// Exponential regression: y = a * e^(b*x). Returns (a, b).
pub fn exponential_regression(x:&[f64],y:&[f64])->(f64,f64){
    let lny:Vec<f64>=y.iter().map(|&yi|crate::exponential::ln(yi)).collect();
    let(ln_a,b)=linear_regression(x,&lny);
    (crate::exponential::exp(ln_a),b)
}

/// Power regression: y = a * x^b. Returns (a, b).
pub fn power_regression(x:&[f64],y:&[f64])->(f64,f64){
    let lnx:Vec<f64>=x.iter().map(|&xi|crate::exponential::ln(xi)).collect();
    let lny:Vec<f64>=y.iter().map(|&yi|crate::exponential::ln(yi)).collect();
    let(ln_a,b)=linear_regression(&lnx,&lny);
    (crate::exponential::exp(ln_a),b)
}

/// Least squares residuals.
pub fn residuals(x:&[f64],y:&[f64])->Vec<f64>{
    let(a,b)=linear_regression(x,y);
    x.iter().zip(y.iter()).map(|(&xi,&yi)|yi-(a+b*xi)).collect()
}

/// Standard error of estimate.
pub fn standard_error(x:&[f64],y:&[f64])->f64{
    let r=residuals(x,y);let n=r.len()as f64;
    powers::sqrt(r.iter().map(|&ri|ri*ri).sum::<f64>()/(n-2.0))
}

/// Moving linear regression (rolling window).
pub fn rolling_regression(x:&[f64],y:&[f64],window:usize)->Vec<(f64,f64)>{
    if x.len()<window{return vec![];}
    (0..=x.len()-window).map(|i|linear_regression(&x[i..i+window],&y[i..i+window])).collect()
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_linear(){let x=[1.0,2.0,3.0,4.0,5.0];let y=[2.0,4.0,6.0,8.0,10.0];
        let(a,b)=linear_regression(&x,&y);assert!(approx_eq(a,0.0,1e-10));assert!(approx_eq(b,2.0,1e-10));}
    #[test]fn test_r_squared(){let x=[1.0,2.0,3.0,4.0,5.0];let y=[2.0,4.0,6.0,8.0,10.0];
        assert!(approx_eq(r_squared(&x,&y),1.0,1e-10));}
    #[test]fn test_poly_regression(){let x=[0.0,1.0,2.0,3.0];let y=[1.0,2.0,5.0,10.0];
        let c=polynomial_regression(&x,&y,2);assert!(approx_eq(poly_predict(&c,0.0),1.0,0.1));}
}
