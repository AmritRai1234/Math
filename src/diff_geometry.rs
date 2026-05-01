//! Differential geometry: Christoffel symbols, geodesics, curvature tensors.

use crate::powers;

/// Christoffel symbols Γ^k_{ij} from metric tensor g.
/// g is n×n, returns n×n×n array.
pub fn christoffel(g:&[Vec<f64>],dg:&Vec<Vec<Vec<f64>>>)->Vec<Vec<Vec<f64>>>{
    let n=g.len();let ginv=crate::matrix::inverse(&g.to_vec()).unwrap();
    let mut gamma=vec![vec![vec![0.0;n];n];n];
    for k in 0..n{for i in 0..n{for j in 0..n{
        for l in 0..n{gamma[k][i][j]+=0.5*ginv[k][l]*(dg[l][i][j]+dg[l][j][i]-dg[i][j][l]);}}}}
    gamma
}

/// Riemann curvature tensor R^l_{ijk} (simplified 2D computation).
pub fn riemann_2d(gamma:&[Vec<Vec<f64>>],dgamma:&[Vec<Vec<Vec<f64>>>])->Vec<Vec<Vec<Vec<f64>>>>{
    let n=gamma.len();
    let mut r=vec![vec![vec![vec![0.0;n];n];n];n];
    for l in 0..n{for i in 0..n{for j in 0..n{for k in 0..n{
        r[l][i][j][k]=dgamma[l][i][k][j]-dgamma[l][i][j][k];
        for m in 0..n{r[l][i][j][k]+=gamma[l][m][j]*gamma[m][i][k]-gamma[l][m][k]*gamma[m][i][j];}}}}}
    r
}

/// Gaussian curvature for 2D surface with metric g.
pub fn gaussian_curvature_2d(r1212:f64,g:&[[f64;2];2])->f64{
    let det=g[0][0]*g[1][1]-g[0][1]*g[1][0];if det.abs()<1e-30{return 0.0;}r1212/det
}

/// Geodesic equation: d²x^k/dt² + Γ^k_{ij} dx^i/dt dx^j/dt = 0.
/// Solve numerically. Returns path as Vec of position vectors.
pub fn geodesic(gamma:&dyn Fn(&[f64])->Vec<Vec<Vec<f64>>>,x0:&[f64],v0:&[f64],dt:f64,steps:usize)->Vec<Vec<f64>>{
    let n=x0.len();let mut x=x0.to_vec();let mut v=v0.to_vec();let mut path=vec![x.clone()];
    for _ in 0..steps{let g=gamma(&x);let mut a=vec![0.0;n];
        for k in 0..n{for i in 0..n{for j in 0..n{a[k]-=g[k][i][j]*v[i]*v[j];}}}
        for i in 0..n{v[i]+=a[i]*dt;x[i]+=v[i]*dt;}path.push(x.clone());}
    path
}

/// First fundamental form coefficients (E, F, G) for parametric surface.
pub fn first_fundamental_form(ru:&[f64;3],rv:&[f64;3])->(f64,f64,f64){
    let e=ru[0]*ru[0]+ru[1]*ru[1]+ru[2]*ru[2];
    let f=ru[0]*rv[0]+ru[1]*rv[1]+ru[2]*rv[2];
    let g=rv[0]*rv[0]+rv[1]*rv[1]+rv[2]*rv[2];(e,f,g)
}

/// Surface area element dA = sqrt(EG - F²) du dv.
pub fn area_element(e:f64,f_coeff:f64,g:f64)->f64{powers::sqrt(e*g-f_coeff*f_coeff)}

/// Mean curvature H = (eG - 2fF + gE) / (2(EG - F²)).
pub fn mean_curvature(e_coeff:f64,f_coeff:f64,g_coeff:f64,l:f64,m:f64,n:f64)->f64{
    let denom=2.0*(e_coeff*g_coeff-f_coeff*f_coeff);
    if denom.abs()<1e-30{return 0.0;}(l*g_coeff-2.0*m*f_coeff+n*e_coeff)/denom
}

/// Principal curvatures from mean curvature H and Gaussian curvature K.
pub fn principal_curvatures(h:f64,k:f64)->(f64,f64){
    let disc=powers::sqrt((h*h-k).max(0.0));(h+disc,h-disc)
}

/// Parallel transport of vector along path (simplified).
pub fn parallel_transport(gamma_fn:&dyn Fn(&[f64])->Vec<Vec<Vec<f64>>>,path:&[Vec<f64>],v0:&[f64])->Vec<f64>{
    let n=v0.len();let mut v=v0.to_vec();
    for step in 0..path.len()-1{let g=gamma_fn(&path[step]);
        let dx:Vec<f64>=(0..n).map(|i|path[step+1][i]-path[step][i]).collect();
        let mut dv=vec![0.0;n];
        for k in 0..n{for i in 0..n{for j in 0..n{dv[k]-=g[k][i][j]*v[i]*dx[j];}}}
        for i in 0..n{v[i]+=dv[i];}}v
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_flat_metric(){// Flat space: Christoffel symbols should be zero
        let g=vec![vec![1.0,0.0],vec![0.0,1.0]];
        let dg=vec![vec![vec![0.0;2];2];2];// dg[l][i][j] = ∂g_{ij}/∂x^l, all zero
        let gamma=christoffel(&g,&dg);
        assert!(approx_eq(gamma[0][0][0],0.0,1e-10));}
    #[test]fn test_fundamental_form(){let ru=[1.0,0.0,0.0];let rv=[0.0,1.0,0.0];
        let(e,f,g)=first_fundamental_form(&ru,&rv);
        assert!(approx_eq(e,1.0,1e-10));assert!(approx_eq(f,0.0,1e-10));assert!(approx_eq(g,1.0,1e-10));}
}
