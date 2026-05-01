//! Linear algebra decompositions: LU, Cholesky, QR, eigenvalues.

use crate::powers;use crate::arithmetic::abs;

/// LU decomposition with partial pivoting. Returns (L, U, permutation).
pub fn lu_decompose(a:&[Vec<f64>])->(Vec<Vec<f64>>,Vec<Vec<f64>>,Vec<usize>){
    let n=a.len();let mut u:Vec<Vec<f64>>=a.to_vec();
    let mut l=vec![vec![0.0;n];n];let mut perm:Vec<usize>=(0..n).collect();
    for k in 0..n{// Partial pivoting
        let mut max_val=abs(u[k][k]);let mut max_row=k;
        for i in k+1..n{if abs(u[i][k])>max_val{max_val=abs(u[i][k]);max_row=i;}}
        if max_row!=k{u.swap(k,max_row);l.swap(k,max_row);perm.swap(k,max_row);}
        l[k][k]=1.0;
        for i in k+1..n{if abs(u[k][k])>1e-30{l[i][k]=u[i][k]/u[k][k];}
            for j in k..n{u[i][j]-=l[i][k]*u[k][j];}}}
    (l,u,perm)
}

/// Solve Ax=b using LU decomposition.
pub fn lu_solve(a:&[Vec<f64>],b:&[f64])->Vec<f64>{
    let n=a.len();let(l,u,perm)=lu_decompose(a);
    let pb:Vec<f64>=perm.iter().map(|&i|b[i]).collect();
    // Forward substitution Ly = Pb
    let mut y=vec![0.0;n];for i in 0..n{y[i]=pb[i];for j in 0..i{y[i]-=l[i][j]*y[j];}}
    // Back substitution Ux = y
    let mut x=vec![0.0;n];for i in (0..n).rev(){x[i]=y[i];for j in i+1..n{x[i]-=u[i][j]*x[j];}
        if abs(u[i][i])>1e-30{x[i]/=u[i][i];}}x
}

/// Cholesky decomposition: A = LL^T (A must be positive definite).
pub fn cholesky(a:&[Vec<f64>])->Option<Vec<Vec<f64>>>{
    let n=a.len();let mut l=vec![vec![0.0;n];n];
    for i in 0..n{for j in 0..=i{let mut sum=0.0;
        for k in 0..j{sum+=l[i][k]*l[j][k];}
        if i==j{let diag=a[i][i]-sum;if diag<=0.0{return None;}l[i][j]=powers::sqrt(diag);}
        else{l[i][j]=(a[i][j]-sum)/l[j][j];}}}Some(l)
}

/// QR decomposition using Gram-Schmidt.
pub fn qr_decompose(a:&[Vec<f64>])->(Vec<Vec<f64>>,Vec<Vec<f64>>){
    let(m,n)=(a.len(),a[0].len());
    let mut q=vec![vec![0.0;m];n]; // columns of Q
    let mut r=vec![vec![0.0;n];n];
    for j in 0..n{// Get column j of A
        let mut v:Vec<f64>=(0..m).map(|i|a[i][j]).collect();
        for i in 0..j{let dot:f64=(0..m).map(|k|q[i][k]*v[k]).sum();r[i][j]=dot;
            for k in 0..m{v[k]-=dot*q[i][k];}}
        let norm=powers::sqrt(v.iter().map(|&x|x*x).sum::<f64>());r[j][j]=norm;
        if norm>1e-30{for k in 0..m{q[j][k]=v[k]/norm;}}}
    // Convert Q from column-major to row-major
    let q_mat:Vec<Vec<f64>>=(0..m).map(|i|(0..n).map(|j|q[j][i]).collect()).collect();
    let r_mat:Vec<Vec<f64>>=r;(q_mat,r_mat)
}

/// Eigenvalues via QR algorithm (real symmetric matrices).
pub fn eigenvalues_qr(a:&[Vec<f64>],max_iter:usize)->Vec<f64>{
    let n=a.len();let mut ak=a.to_vec();
    for _ in 0..max_iter{let(q,r)=qr_decompose(&ak);
        // A_{k+1} = R * Q
        ak=vec![vec![0.0;n];n];
        for i in 0..n{for j in 0..n{for k in 0..n{ak[i][j]+=r[i][k]*q[k][j];}}}}
    (0..n).map(|i|ak[i][i]).collect()
}

/// Solve Ax=b using Cholesky (A positive definite).
pub fn cholesky_solve(a:&[Vec<f64>],b:&[f64])->Option<Vec<f64>>{
    let l=cholesky(a)?;let n=a.len();
    // Forward: Ly = b
    let mut y=vec![0.0;n];for i in 0..n{y[i]=b[i];for j in 0..i{y[i]-=l[i][j]*y[j];}y[i]/=l[i][i];}
    // Backward: L^T x = y
    let mut x=vec![0.0;n];for i in (0..n).rev(){x[i]=y[i];for j in i+1..n{x[i]-=l[j][i]*x[j];}x[i]/=l[i][i];}
    Some(x)
}

/// Condition number estimate (ratio of max/min singular values via eigenvalues of A^T A).
pub fn condition_number(a:&[Vec<f64>])->f64{
    let n=a.len();let m=a[0].len();
    let mut ata=vec![vec![0.0;m];m];
    for i in 0..m{for j in 0..m{for k in 0..n{ata[i][j]+=a[k][i]*a[k][j];}}}
    let eigs=eigenvalues_qr(&ata,100);
    let max=eigs.iter().cloned().fold(0.0f64,f64::max);
    let min=eigs.iter().cloned().fold(f64::INFINITY,f64::min);
    if min.abs()<1e-30{f64::INFINITY}else{powers::sqrt(max/min)}
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_lu_solve(){let a=vec![vec![2.0,1.0],vec![1.0,3.0]];let b=vec![5.0,7.0];
        let x=lu_solve(&a,&b);assert!(approx_eq(2.0*x[0]+x[1],5.0,1e-8));}
    #[test]fn test_cholesky(){let a=vec![vec![4.0,2.0],vec![2.0,3.0]];let l=cholesky(&a).unwrap();
        assert!(approx_eq(l[0][0],2.0,1e-10));}
    #[test]fn test_qr(){let a=vec![vec![1.0,1.0],vec![0.0,1.0],vec![1.0,0.0]];
        let(q,r)=qr_decompose(&a);assert!(approx_eq(r[0][0],powers::sqrt(2.0),1e-8));}
}
