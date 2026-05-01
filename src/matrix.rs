//! Advanced matrix operations: QR, eigenvalues, SVD, rank, norms.

use crate::arithmetic::abs;
use crate::powers;

type Matrix = Vec<Vec<f64>>;

/// Matrix-vector multiply.
pub fn mat_vec_mul(a:&Matrix,v:&[f64])->Vec<f64>{
    a.iter().map(|row|row.iter().zip(v.iter()).map(|(&a,&b)|a*b).sum()).collect()
}

/// Matrix-matrix multiply.
pub fn mat_mul(a:&Matrix,b:&Matrix)->Matrix{
    let m=a.len();let p=b[0].len();let n=b.len();
    (0..m).map(|i|(0..p).map(|j|(0..n).map(|k|a[i][k]*b[k][j]).sum()).collect()).collect()
}

/// Transpose.
pub fn transpose(a:&Matrix)->Matrix{
    let m=a.len();let n=a[0].len();
    (0..n).map(|j|(0..m).map(|i|a[i][j]).collect()).collect()
}

/// Frobenius norm.
pub fn frobenius_norm(a:&Matrix)->f64{
    powers::sqrt(a.iter().flat_map(|r|r.iter()).map(|&x|x*x).sum())
}

/// Infinity norm (max row sum).
pub fn inf_norm(a:&Matrix)->f64{
    a.iter().map(|r|r.iter().map(|&x|abs(x)).sum::<f64>()).fold(0.0f64,|a,b|a.max(b))
}

/// Matrix rank via row echelon form.
pub fn rank(a:&Matrix)->usize{
    let mut m=a.clone();let rows=m.len();if rows==0{return 0;}let cols=m[0].len();
    let mut r=0;
    for c in 0..cols{
        if r>=rows{break;}
        let mut pivot=None;for i in r..rows{if abs(m[i][c])>1e-12{pivot=Some(i);break;}}
        if let Some(p)=pivot{m.swap(r,p);let scale=m[r][c];
            for j in 0..cols{m[r][j]/=scale;}
            for i in 0..rows{if i!=r{let f=m[i][c];for j in 0..cols{m[i][j]-=f*m[r][j];}}}
            r+=1;
        }
    }
    r
}

/// QR decomposition via Gram-Schmidt. Returns (Q, R).
pub fn qr_decompose(a:&Matrix)->(Matrix,Matrix){
    let m=a.len();let n=a[0].len();let k=m.min(n);
    let mut q:Vec<Vec<f64>>=vec![vec![0.0;m];k];
    let mut r=vec![vec![0.0;n];k];

    for j in 0..k{
        // Get column j of A
        let mut v:Vec<f64>=(0..m).map(|i|a[i][j]).collect();
        // Subtract projections
        for i in 0..j{let dot:f64=(0..m).map(|row|q[i][row]*v[row]).sum();
            r[i][j]=dot;for row in 0..m{v[row]-=dot*q[i][row];}}
        let norm=powers::sqrt(v.iter().map(|&x|x*x).sum());
        r[j][j]=norm;
        if norm>1e-15{for row in 0..m{q[j][row]=v[row]/norm;}}
    }

    // Convert Q from row-major basis vectors to column-major matrix
    let q_mat:Matrix=(0..m).map(|i|(0..k).map(|j|q[j][i]).collect()).collect();
    let r_mat:Matrix=(0..k).map(|i|(0..n).map(|j|r[i][j]).collect()).collect();
    (q_mat,r_mat)
}

/// Eigenvalues of a symmetric matrix using QR algorithm.
pub fn eigenvalues_symmetric(a:&Matrix)->Vec<f64>{
    let n=a.len();let mut m=a.clone();
    for _ in 0..2000{let(q,r)=qr_decompose(&m);
        // A_{k+1} = R * Q (not R * Q^T)
        let q_square:Matrix=(0..n).map(|i|(0..n).map(|j|q[i][j]).collect()).collect();
        m=mat_mul(&r,&q_square);
        let mut off=0.0;for i in 0..n{for j in 0..n{if i!=j{off+=m[i][j]*m[i][j];}}}
        if off<1e-20{break;}}
    (0..n).map(|i|m[i][i]).collect()
}

/// Cholesky decomposition of positive-definite matrix. Returns L where A = L*L^T.
pub fn cholesky(a:&Matrix)->Option<Matrix>{
    let n=a.len();let mut l=vec![vec![0.0;n];n];
    for i in 0..n{for j in 0..=i{
        let mut sum=0.0;for k in 0..j{sum+=l[i][k]*l[j][k];}
        if i==j{let val=a[i][i]-sum;if val<=0.0{return None;}l[i][j]=powers::sqrt(val);}
        else{l[i][j]=(a[i][j]-sum)/l[j][j];}
    }}
    Some(l)
}

/// Determinant via LU decomposition (for any size).
pub fn determinant(a:&Matrix)->f64{
    let n=a.len();let mut m=a.clone();let mut det=1.0;
    for i in 0..n{
        let mut pivot=i;for j in(i+1)..n{if abs(m[j][i])>abs(m[pivot][i]){pivot=j;}}
        if abs(m[pivot][i])<1e-15{return 0.0;}
        if pivot!=i{m.swap(i,pivot);det=-det;}
        det*=m[i][i];
        for j in(i+1)..n{let f=m[j][i]/m[i][i];for k in i..n{m[j][k]-=f*m[i][k];}}
    }
    det
}

/// Matrix inverse via Gauss-Jordan.
pub fn inverse(a:&Matrix)->Option<Matrix>{
    let n=a.len();let mut aug:Matrix=a.iter().enumerate().map(|(i,row)|{
        let mut r=row.clone();r.resize(2*n,0.0);r[n+i]=1.0;r}).collect();
    for i in 0..n{
        let mut pivot=i;for j in(i+1)..n{if abs(aug[j][i])>abs(aug[pivot][i]){pivot=j;}}
        if abs(aug[pivot][i])<1e-15{return None;}
        aug.swap(i,pivot);let scale=aug[i][i];
        for j in 0..2*n{aug[i][j]/=scale;}
        for j in 0..n{if j!=i{let f=aug[j][i];for k in 0..2*n{aug[j][k]-=f*aug[i][k];}}}
    }
    Some(aug.iter().map(|r|r[n..].to_vec()).collect())
}

/// Condition number (ratio of largest to smallest singular value, approximated).
pub fn condition_number(a:&Matrix)->f64{
    let ata=mat_mul(&transpose(a),a);
    let eigs=eigenvalues_symmetric(&ata);
    let max_eig=eigs.iter().cloned().fold(0.0f64,f64::max);
    let min_eig=eigs.iter().cloned().fold(f64::MAX,f64::min);
    if min_eig<=0.0{return f64::INFINITY;}
    powers::sqrt(max_eig/min_eig)
}

/// Kronecker product.
pub fn kronecker(a:&Matrix,b:&Matrix)->Matrix{
    let(ma,na)=(a.len(),a[0].len());let(mb,nb)=(b.len(),b[0].len());
    let mut result=vec![vec![0.0;na*nb];ma*mb];
    for i in 0..ma{for j in 0..na{for k in 0..mb{for l in 0..nb{
        result[i*mb+k][j*nb+l]=a[i][j]*b[k][l];}}}}
    result
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_qr(){
        let a=vec![vec![1.0,1.0],vec![0.0,1.0],vec![1.0,0.0]];
        let(q,r)=qr_decompose(&a);let reconstructed=mat_mul(&q,&r);
        for i in 0..3{for j in 0..2{assert!(approx_eq(reconstructed[i][j],a[i][j],1e-10));}}
    }
    #[test]fn test_cholesky(){
        let a=vec![vec![4.0,2.0],vec![2.0,3.0]];let l=cholesky(&a).unwrap();
        let lt=transpose(&l);let result=mat_mul(&l,&lt);
        assert!(approx_eq(result[0][0],4.0,1e-10));assert!(approx_eq(result[0][1],2.0,1e-10));
    }
    #[test]fn test_det(){
        let a=vec![vec![1.0,2.0],vec![3.0,4.0]];assert!(approx_eq(determinant(&a),-2.0,1e-10));
    }
    #[test]fn test_inverse(){
        let a=vec![vec![1.0,2.0],vec![3.0,4.0]];let inv=inverse(&a).unwrap();
        let id=mat_mul(&a,&inv);assert!(approx_eq(id[0][0],1.0,1e-10));assert!(approx_eq(id[0][1],0.0,1e-10));
    }
    #[test]fn test_rank(){
        let a=vec![vec![1.0,2.0,3.0],vec![4.0,5.0,6.0],vec![7.0,8.0,9.0]];
        assert_eq!(rank(&a),2);// rank-deficient
    }
    #[test]fn test_eigenvalues(){
        let a=vec![vec![2.0,1.0],vec![1.0,2.0]];
        let mut eigs=eigenvalues_symmetric(&a);eigs.sort_by(|a,b|a.partial_cmp(b).unwrap());
        assert!(approx_eq(eigs[0],1.0,1e-6));assert!(approx_eq(eigs[1],3.0,1e-6));
    }
}
