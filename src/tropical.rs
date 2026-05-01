//! Tropical algebra (min-plus and max-plus semirings).

/// Tropical addition (min): a ⊕ b = min(a, b).
pub fn trop_add_min(a:f64,b:f64)->f64{a.min(b)}
/// Tropical multiplication: a ⊗ b = a + b.
pub fn trop_mul(a:f64,b:f64)->f64{a+b}
/// Tropical zero (identity for ⊕ in min-plus): ∞.
pub const TROP_ZERO:f64=f64::INFINITY;
/// Tropical one (identity for ⊗): 0.
pub const TROP_ONE:f64=0.0;

/// Max-plus tropical addition: a ⊕ b = max(a, b).
pub fn trop_add_max(a:f64,b:f64)->f64{a.max(b)}

/// Tropical matrix multiplication (min-plus).
pub fn trop_mat_mul(a:&[Vec<f64>],b:&[Vec<f64>])->Vec<Vec<f64>>{
    let m=a.len();let n=b[0].len();let p=b.len();
    (0..m).map(|i|(0..n).map(|j|(0..p).map(|k|trop_mul(a[i][k],b[k][j])).fold(TROP_ZERO,trop_add_min)).collect()).collect()
}

/// Tropical matrix power (for shortest paths).
pub fn trop_mat_pow(a:&[Vec<f64>],exp:u32)->Vec<Vec<f64>>{
    let n=a.len();let mut result=vec![vec![TROP_ZERO;n];n];
    for i in 0..n{result[i][i]=TROP_ONE;}
    let mut base=a.to_vec();let mut e=exp;
    while e>0{if e&1==1{result=trop_mat_mul(&result,&base);}base=trop_mat_mul(&base,&base);e>>=1;}result
}

/// Tropical convex hull (min of linear functions).
pub fn trop_convex_eval(coeffs:&[(f64,f64)],x:f64)->f64{
    coeffs.iter().map(|&(a,b)|trop_mul(a,x)+b).fold(TROP_ZERO,trop_add_min)
}

/// Tropical polynomial evaluation: min(a_i + i*x).
pub fn trop_poly_eval(coeffs:&[f64],x:f64)->f64{
    coeffs.iter().enumerate().map(|(i,&c)|c+i as f64*x).fold(TROP_ZERO,trop_add_min)
}

/// Tropical determinant (min of sum over all permutations).
pub fn trop_det(a:&[Vec<f64>])->f64{
    let n=a.len();if n==1{return a[0][0];}
    let perms=crate::abstract_algebra::all_permutations(n);
    perms.iter().map(|p|{(0..n).map(|i|a[i][p.apply(i)]).sum::<f64>()}).fold(TROP_ZERO,trop_add_min)
}

/// Kleene star (transitive closure) for tropical matrices.
pub fn trop_kleene_star(a:&[Vec<f64>])->Vec<Vec<f64>>{
    let n=a.len();let mut result=vec![vec![TROP_ZERO;n];n];
    for i in 0..n{result[i][i]=TROP_ONE;}
    for k in 0..n{for i in 0..n{for j in 0..n{
        let through_k=trop_mul(result[i][k],trop_mul(a[k][k].min(TROP_ONE),result[k][j]));
        result[i][j]=trop_add_min(result[i][j],through_k);}}}
    result
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_add(){assert!(approx_eq(trop_add_min(3.0,5.0),3.0,1e-10));}
    #[test]fn test_mul(){assert!(approx_eq(trop_mul(3.0,5.0),8.0,1e-10));}
    #[test]fn test_identity(){assert!(approx_eq(trop_add_min(3.0,TROP_ZERO),3.0,1e-10));assert!(approx_eq(trop_mul(3.0,TROP_ONE),3.0,1e-10));}
    #[test]fn test_mat_mul(){let a=vec![vec![0.0,1.0],vec![2.0,TROP_ZERO]];let b=vec![vec![3.0,TROP_ZERO],vec![TROP_ZERO,4.0]];
        let c=trop_mat_mul(&a,&b);assert!(approx_eq(c[0][0],3.0,1e-10));assert!(approx_eq(c[0][1],5.0,1e-10));}
}
