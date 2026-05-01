//! Basic tensor operations.

use crate::arithmetic::abs;

/// 3D tensor (rank 3).
#[derive(Clone,Debug)]
pub struct Tensor3{pub data:Vec<Vec<Vec<f64>>>,pub shape:(usize,usize,usize)}

impl Tensor3{
    pub fn zeros(d0:usize,d1:usize,d2:usize)->Self{Self{data:vec![vec![vec![0.0;d2];d1];d0],shape:(d0,d1,d2)}}
    pub fn get(&self,i:usize,j:usize,k:usize)->f64{self.data[i][j][k]}
    pub fn set(&mut self,i:usize,j:usize,k:usize,v:f64){self.data[i][j][k]=v;}
    pub fn add(&self,o:&Self)->Self{let mut r=self.clone();for i in 0..self.shape.0{for j in 0..self.shape.1{for k in 0..self.shape.2{r.data[i][j][k]+=o.data[i][j][k];}}}r}
    pub fn scale(&self,s:f64)->Self{let mut r=self.clone();for i in 0..self.shape.0{for j in 0..self.shape.1{for k in 0..self.shape.2{r.data[i][j][k]*=s;}}}r}
    pub fn norm(&self)->f64{crate::powers::sqrt(self.data.iter().flat_map(|a|a.iter().flat_map(|b|b.iter())).map(|&x|x*x).sum())}
    /// Contract indices j and k (trace over last two indices).
    pub fn contract_jk(&self)->Vec<f64>{
        let n=self.shape.1.min(self.shape.2);
        (0..self.shape.0).map(|i|(0..n).map(|j|self.data[i][j][j]).sum()).collect()
    }
}

/// Levi-Civita symbol ε_{ijk} for 3D.
pub fn levi_civita(i:usize,j:usize,k:usize)->i8{
    if i==j||j==k||i==k{return 0;}
    // Even permutations of (0,1,2): (0,1,2),(1,2,0),(2,0,1) → +1
    // Odd permutations: (0,2,1),(1,0,2),(2,1,0) → -1
    match(i,j,k){(0,1,2)|(1,2,0)|(2,0,1)=>1,(0,2,1)|(2,1,0)|(1,0,2)=>-1,_=>0}
}

/// Kronecker delta δ_{ij}.
pub fn kronecker_delta(i:usize,j:usize)->f64{if i==j{1.0}else{0.0}}

/// Tensor outer product of two vectors: T_{ij} = a_i * b_j.
pub fn outer_product_2(a:&[f64],b:&[f64])->Vec<Vec<f64>>{
    a.iter().map(|&ai|b.iter().map(|&bj|ai*bj).collect()).collect()
}

/// Tensor outer product of three vectors: T_{ijk} = a_i * b_j * c_k.
pub fn outer_product_3(a:&[f64],b:&[f64],c:&[f64])->Tensor3{
    let mut t=Tensor3::zeros(a.len(),b.len(),c.len());
    for i in 0..a.len(){for j in 0..b.len(){for k in 0..c.len(){t.data[i][j][k]=a[i]*b[j]*c[k];}}}t
}

/// Einstein summation: C_{ik} = Σ_j A_{ij} * B_{jk} (matrix multiply as tensor contraction).
pub fn einsum_ij_jk(a:&[Vec<f64>],b:&[Vec<f64>])->Vec<Vec<f64>>{crate::matrix::mat_mul(&a.to_vec(),&b.to_vec())}

/// Tensor contraction: contract a rank-3 tensor over index 1 with a vector.
pub fn contract_with_vector(t:&Tensor3,v:&[f64])->Vec<Vec<f64>>{
    let(d0,d1,d2)=t.shape;
    (0..d0).map(|i|(0..d2).map(|k|(0..d1).map(|j|t.data[i][j][k]*v[j]).sum()).collect()).collect()
}

/// Metric tensor for spherical coordinates (r, θ, φ).
pub fn spherical_metric(r:f64,theta:f64)->[[f64;3];3]{
    let st=crate::trig::sin(theta);
    [[1.0,0.0,0.0],[0.0,r*r,0.0],[0.0,0.0,r*r*st*st]]
}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_levi_civita(){assert_eq!(levi_civita(0,1,2),1);assert_eq!(levi_civita(1,0,2),-1);assert_eq!(levi_civita(0,0,1),0);}
    #[test]fn test_kronecker(){assert_eq!(kronecker_delta(1,1),1.0);assert_eq!(kronecker_delta(0,1),0.0);}
    #[test]fn test_outer(){let a=[1.0,2.0];let b=[3.0,4.0];let t=outer_product_2(&a,&b);assert_eq!(t[0][0],3.0);assert_eq!(t[1][1],8.0);}
    #[test]fn test_tensor3(){let t=Tensor3::zeros(2,2,2);assert_eq!(t.norm(),0.0);}
}
