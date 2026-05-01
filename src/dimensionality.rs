//! Dimensionality reduction: PCA, SVD-based methods.

use crate::matrix;
use crate::stats;
use crate::powers;

/// Principal Component Analysis. Returns (components, explained_variance, projections).
pub fn pca(data:&[Vec<f64>],n_components:usize)->(Vec<Vec<f64>>,Vec<f64>,Vec<Vec<f64>>){
    let n=data.len();let d=data[0].len();
    // Center data
    let means:Vec<f64>=(0..d).map(|j|data.iter().map(|row|row[j]).sum::<f64>()/n as f64).collect();
    let centered:Vec<Vec<f64>>=data.iter().map(|row|row.iter().zip(means.iter()).map(|(&x,&m)|x-m).collect()).collect();
    // Covariance matrix
    let mut cov=vec![vec![0.0;d];d];
    for i in 0..d{for j in 0..d{cov[i][j]=centered.iter().map(|row|row[i]*row[j]).sum::<f64>()/(n-1)as f64;}}
    // Eigendecomposition
    let eigenvalues=matrix::eigenvalues_symmetric(&cov);
    // Get eigenvectors via QR
    let(q,_)=matrix::qr_decompose(&cov);
    // Sort by eigenvalue (descending)
    let mut eig_pairs:Vec<(f64,Vec<f64>)>=eigenvalues.iter().enumerate().map(|(i,&val)|{
        let vec:Vec<f64>=(0..q.len()).map(|r|q[r][i]).collect();(val,vec)}).collect();
    eig_pairs.sort_by(|a,b|b.0.partial_cmp(&a.0).unwrap());
    let components:Vec<Vec<f64>>=eig_pairs.iter().take(n_components).map(|(_,v)|v.clone()).collect();
    let explained:Vec<f64>=eig_pairs.iter().take(n_components).map(|(v,_)|*v).collect();
    // Project data
    let projections:Vec<Vec<f64>>=centered.iter().map(|row|{
        components.iter().map(|comp|comp.iter().zip(row.iter()).map(|(&c,&x)|c*x).sum()).collect()}).collect();
    (components,explained,projections)
}

/// Explained variance ratio.
pub fn explained_variance_ratio(eigenvalues:&[f64])->Vec<f64>{
    let total:f64=eigenvalues.iter().sum();eigenvalues.iter().map(|&v|v/total).collect()
}

/// Cumulative explained variance.
pub fn cumulative_variance(eigenvalues:&[f64])->Vec<f64>{
    let ratios=explained_variance_ratio(eigenvalues);
    let mut cum=Vec::new();let mut sum=0.0;
    for r in ratios{sum+=r;cum.push(sum);}cum
}

/// Number of components needed to explain target_variance (0 to 1).
pub fn n_components_for_variance(eigenvalues:&[f64],target:f64)->usize{
    let cum=cumulative_variance(eigenvalues);
    cum.iter().position(|&v|v>=target).map(|p|p+1).unwrap_or(eigenvalues.len())
}

/// Standardize data (zero mean, unit variance).
pub fn standardize(data:&[Vec<f64>])->Vec<Vec<f64>>{
    let n=data.len();let d=data[0].len();
    let means:Vec<f64>=(0..d).map(|j|data.iter().map(|row|row[j]).sum::<f64>()/n as f64).collect();
    let stds:Vec<f64>=(0..d).map(|j|{let m=means[j];
        powers::sqrt(data.iter().map(|row|(row[j]-m)*(row[j]-m)).sum::<f64>()/(n-1)as f64)}).collect();
    data.iter().map(|row|row.iter().enumerate().map(|(j,&x)|if stds[j]>1e-30{(x-means[j])/stds[j]}else{0.0}).collect()).collect()
}

/// Covariance matrix.
pub fn covariance_matrix(data:&[Vec<f64>])->Vec<Vec<f64>>{
    let n=data.len();let d=data[0].len();
    let means:Vec<f64>=(0..d).map(|j|data.iter().map(|row|row[j]).sum::<f64>()/n as f64).collect();
    let mut cov=vec![vec![0.0;d];d];
    for i in 0..d{for j in 0..d{cov[i][j]=data.iter().map(|row|(row[i]-means[i])*(row[j]-means[j])).sum::<f64>()/(n-1)as f64;}}cov
}

/// Correlation matrix.
pub fn correlation_matrix(data:&[Vec<f64>])->Vec<Vec<f64>>{
    let cov=covariance_matrix(data);let d=cov.len();
    let stds:Vec<f64>=(0..d).map(|i|powers::sqrt(cov[i][i])).collect();
    (0..d).map(|i|(0..d).map(|j|if stds[i]>1e-30&&stds[j]>1e-30{cov[i][j]/(stds[i]*stds[j])}else{0.0}).collect()).collect()
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_pca(){let data=vec![vec![1.0,2.0],vec![3.0,4.0],vec![5.0,6.0],vec![7.0,8.0]];
        let(comp,_,proj)=pca(&data,1);assert_eq!(comp.len(),1);assert_eq!(proj.len(),4);assert_eq!(proj[0].len(),1);}
    #[test]fn test_standardize(){let data=vec![vec![1.0,2.0],vec![3.0,4.0],vec![5.0,6.0]];
        let std=standardize(&data);let mean:f64=std.iter().map(|r|r[0]).sum::<f64>()/3.0;assert!(approx_eq(mean,0.0,1e-10));}
    #[test]fn test_corr(){let data=vec![vec![1.0,2.0],vec![2.0,4.0],vec![3.0,6.0]];
        let c=correlation_matrix(&data);assert!(approx_eq(c[0][1],1.0,1e-6));}// perfectly correlated
}
