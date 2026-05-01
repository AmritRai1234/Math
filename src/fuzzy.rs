//! Fuzzy logic: fuzzy sets, membership functions, operations.

use crate::arithmetic;

/// Membership functions.
pub fn triangular(x:f64,a:f64,b:f64,c:f64)->f64{
    if x<=a||x>=c{0.0}else if x<=b{(x-a)/(b-a)}else{(c-x)/(c-b)}
}

pub fn trapezoidal(x:f64,a:f64,b:f64,c:f64,d:f64)->f64{
    if x<=a||x>=d{0.0}else if x>=b&&x<=c{1.0}else if x<b{(x-a)/(b-a)}else{(d-x)/(d-c)}
}

pub fn gaussian(x:f64,mean:f64,sigma:f64)->f64{
    crate::exponential::exp(-0.5*((x-mean)/sigma)*((x-mean)/sigma))
}

pub fn bell(x:f64,a:f64,b:f64,c:f64)->f64{
    1.0/(1.0+crate::powers::pow(arithmetic::abs((x-c)/a),2.0*b))
}

pub fn sigmoid(x:f64,a:f64,c:f64)->f64{
    1.0/(1.0+crate::exponential::exp(-a*(x-c)))
}

/// Fuzzy set: vector of (value, membership degree) pairs.
pub type FuzzySet=Vec<(f64,f64)>;

/// Fuzzy union: max of membership degrees.
pub fn fuzzy_union(a:&FuzzySet,b:&FuzzySet)->FuzzySet{
    a.iter().zip(b.iter()).map(|(&(x,ma),&(_,mb))|(x,ma.max(mb))).collect()
}

/// Fuzzy intersection: min of membership degrees.
pub fn fuzzy_intersection(a:&FuzzySet,b:&FuzzySet)->FuzzySet{
    a.iter().zip(b.iter()).map(|(&(x,ma),&(_,mb))|(x,ma.min(mb))).collect()
}

/// Fuzzy complement.
pub fn fuzzy_complement(a:&FuzzySet)->FuzzySet{
    a.iter().map(|&(x,m)|(x,1.0-m)).collect()
}

/// Alpha-cut: elements with membership ≥ alpha.
pub fn alpha_cut(a:&FuzzySet,alpha:f64)->Vec<f64>{
    a.iter().filter(|&&(_,m)|m>=alpha).map(|&(x,_)|x).collect()
}

/// Fuzzy AND (t-norm: minimum).
pub fn t_norm_min(a:f64,b:f64)->f64{a.min(b)}
/// Fuzzy AND (t-norm: product).
pub fn t_norm_product(a:f64,b:f64)->f64{a*b}
/// Fuzzy AND (t-norm: Łukasiewicz).
pub fn t_norm_luk(a:f64,b:f64)->f64{(a+b-1.0).max(0.0)}

/// Fuzzy OR (s-norm: maximum).
pub fn s_norm_max(a:f64,b:f64)->f64{a.max(b)}
/// Fuzzy OR (s-norm: probabilistic sum).
pub fn s_norm_prob(a:f64,b:f64)->f64{a+b-a*b}

/// Defuzzification: centroid method.
pub fn defuzzify_centroid(set:&FuzzySet)->f64{
    let num:f64=set.iter().map(|&(x,m)|x*m).sum();
    let den:f64=set.iter().map(|&(_,m)|m).sum();
    if den==0.0{0.0}else{num/den}
}

/// Defuzzification: mean of maximum.
pub fn defuzzify_mom(set:&FuzzySet)->f64{
    let max_m=set.iter().map(|&(_,m)|m).fold(0.0f64,f64::max);
    let maxes:Vec<f64>=set.iter().filter(|&&(_,m)|(m-max_m).abs()<1e-10).map(|&(x,_)|x).collect();
    if maxes.is_empty(){0.0}else{maxes.iter().sum::<f64>()/maxes.len()as f64}
}

/// Fuzzy entropy.
pub fn fuzzy_entropy(set:&FuzzySet)->f64{
    let n=set.len()as f64;if n==0.0{return 0.0;}
    set.iter().map(|&(_,m)|{let mn=m.min(1.0-m);let mx=m.max(1.0-m);if mx>0.0{mn/mx}else{0.0}}).sum::<f64>()/n
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_triangular(){assert!(approx_eq(triangular(5.0,0.0,5.0,10.0),1.0,1e-10));assert!(approx_eq(triangular(0.0,0.0,5.0,10.0),0.0,1e-10));}
    #[test]fn test_gaussian_peak(){assert!(approx_eq(gaussian(5.0,5.0,1.0),1.0,1e-10));}
    #[test]fn test_defuzzify(){let set=vec![(1.0,0.2),(2.0,0.5),(3.0,0.8),(4.0,0.5),(5.0,0.2)];
        let c=defuzzify_centroid(&set);assert!(c>2.5&&c<3.5);}
    #[test]fn test_complement(){let s=vec![(1.0,0.3),(2.0,0.7)];let c=fuzzy_complement(&s);
        assert!(approx_eq(c[0].1,0.7,1e-10));assert!(approx_eq(c[1].1,0.3,1e-10));}
}
