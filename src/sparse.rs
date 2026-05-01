//! Sparse matrices in CSR (Compressed Sparse Row) format.

use crate::arithmetic::abs;

#[derive(Clone,Debug)]
pub struct SparseMatrix{pub rows:usize,pub cols:usize,pub row_ptr:Vec<usize>,pub col_idx:Vec<usize>,pub values:Vec<f64>}

impl SparseMatrix{
    pub fn new(rows:usize,cols:usize)->Self{Self{rows,cols,row_ptr:vec![0;rows+1],col_idx:Vec::new(),values:Vec::new()}}

    /// Build from triplets (row, col, value).
    pub fn from_triplets(rows:usize,cols:usize,triplets:&[(usize,usize,f64)])->Self{
        let mut sorted=triplets.to_vec();sorted.sort_by_key(|&(r,c,_)|(r,c));
        let mut row_ptr=vec![0;rows+1];let mut col_idx=Vec::new();let mut values=Vec::new();
        for &(r,c,v) in &sorted{if abs(v)>1e-30{col_idx.push(c);values.push(v);row_ptr[r+1]+=1;}}
        for i in 0..rows{row_ptr[i+1]+=row_ptr[i];}
        Self{rows,cols,row_ptr,col_idx,values}
    }

    /// Build from dense matrix.
    pub fn from_dense(dense:&[Vec<f64>])->Self{
        let rows=dense.len();let cols=if rows>0{dense[0].len()}else{0};
        let mut triplets=Vec::new();
        for(i,row) in dense.iter().enumerate(){for(j,&v) in row.iter().enumerate(){if abs(v)>1e-30{triplets.push((i,j,v));}}}
        Self::from_triplets(rows,cols,&triplets)
    }

    /// Number of non-zero elements.
    pub fn nnz(&self)->usize{self.values.len()}

    /// Density (ratio of non-zeros to total elements).
    pub fn density(&self)->f64{self.nnz()as f64/(self.rows*self.cols)as f64}

    /// Get element at (i, j).
    pub fn get(&self,i:usize,j:usize)->f64{
        for k in self.row_ptr[i]..self.row_ptr[i+1]{if self.col_idx[k]==j{return self.values[k];}}0.0
    }

    /// Sparse matrix-vector multiply: y = A*x.
    pub fn mul_vec(&self,x:&[f64])->Vec<f64>{
        let mut y=vec![0.0;self.rows];
        for i in 0..self.rows{for k in self.row_ptr[i]..self.row_ptr[i+1]{y[i]+=self.values[k]*x[self.col_idx[k]];}}y
    }

    /// Sparse matrix-matrix add.
    pub fn add(&self,other:&Self)->Self{
        let mut triplets=Vec::new();
        for i in 0..self.rows{for k in self.row_ptr[i]..self.row_ptr[i+1]{triplets.push((i,self.col_idx[k],self.values[k]));}
            for k in other.row_ptr[i]..other.row_ptr[i+1]{triplets.push((i,other.col_idx[k],other.values[k]));}}
        // Merge duplicates
        triplets.sort_by_key(|&(r,c,_)|(r,c));
        let mut merged=Vec::new();let mut i=0;
        while i<triplets.len(){let(r,c,mut v)=triplets[i];i+=1;
            while i<triplets.len()&&triplets[i].0==r&&triplets[i].1==c{v+=triplets[i].2;i+=1;}
            if abs(v)>1e-30{merged.push((r,c,v));}}
        Self::from_triplets(self.rows,self.cols,&merged)
    }

    /// Scale by constant.
    pub fn scale(&self,s:f64)->Self{let v:Vec<f64>=self.values.iter().map(|&x|x*s).collect();
        Self{rows:self.rows,cols:self.cols,row_ptr:self.row_ptr.clone(),col_idx:self.col_idx.clone(),values:v}}

    /// Convert to dense.
    pub fn to_dense(&self)->Vec<Vec<f64>>{
        let mut d=vec![vec![0.0;self.cols];self.rows];
        for i in 0..self.rows{for k in self.row_ptr[i]..self.row_ptr[i+1]{d[i][self.col_idx[k]]=self.values[k];}}d
    }

    /// Transpose.
    pub fn transpose(&self)->Self{
        let mut triplets=Vec::new();
        for i in 0..self.rows{for k in self.row_ptr[i]..self.row_ptr[i+1]{triplets.push((self.col_idx[k],i,self.values[k]));}}
        Self::from_triplets(self.cols,self.rows,&triplets)
    }

    /// Sparse matrix-matrix multiply.
    pub fn mul_mat(&self,other:&Self)->Self{
        let mut triplets=Vec::new();
        for i in 0..self.rows{for k in self.row_ptr[i]..self.row_ptr[i+1]{let j=self.col_idx[k];let a=self.values[k];
            for l in other.row_ptr[j]..other.row_ptr[j+1]{triplets.push((i,other.col_idx[l],a*other.values[l]));}}}
        let mut sorted=triplets;sorted.sort_by_key(|&(r,c,_)|(r,c));
        let mut merged=Vec::new();let mut idx=0;
        while idx<sorted.len(){let(r,c,mut v)=sorted[idx];idx+=1;
            while idx<sorted.len()&&sorted[idx].0==r&&sorted[idx].1==c{v+=sorted[idx].2;idx+=1;}
            if abs(v)>1e-30{merged.push((r,c,v));}}
        Self::from_triplets(self.rows,other.cols,&merged)
    }

    /// Conjugate gradient solver: Ax = b for symmetric positive-definite A.
    pub fn solve_cg(&self,b:&[f64],tol:f64,max_iter:usize)->Vec<f64>{
        let n=self.rows;let mut x=vec![0.0;n];
        let mut r=b.to_vec();let mut p=r.clone();let mut rs:f64=r.iter().map(|&x|x*x).sum();
        for _ in 0..max_iter{let ap=self.mul_vec(&p);
            let pap:f64=p.iter().zip(ap.iter()).map(|(&a,&b)|a*b).sum();
            if abs(pap)<1e-30{break;}let alpha=rs/pap;
            for i in 0..n{x[i]+=alpha*p[i];r[i]-=alpha*ap[i];}
            let rs_new:f64=r.iter().map(|&x|x*x).sum();
            if rs_new<tol*tol{break;}let beta=rs_new/rs;
            for i in 0..n{p[i]=r[i]+beta*p[i];}rs=rs_new;}
        x
    }

    /// Diagonal elements.
    pub fn diagonal(&self)->Vec<f64>{(0..self.rows.min(self.cols)).map(|i|self.get(i,i)).collect()}

    /// Frobenius norm.
    pub fn frobenius_norm(&self)->f64{crate::powers::sqrt(self.values.iter().map(|&x|x*x).sum())}
}

/// Create sparse identity matrix.
pub fn sparse_identity(n:usize)->SparseMatrix{
    let triplets:Vec<(usize,usize,f64)>=(0..n).map(|i|(i,i,1.0)).collect();
    SparseMatrix::from_triplets(n,n,&triplets)
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_mul_vec(){let m=SparseMatrix::from_triplets(2,2,&[(0,0,2.0),(0,1,1.0),(1,0,1.0),(1,1,3.0)]);
        let y=m.mul_vec(&[1.0,2.0]);assert!(approx_eq(y[0],4.0,1e-10));assert!(approx_eq(y[1],7.0,1e-10));}
    #[test]fn test_transpose(){let m=SparseMatrix::from_triplets(2,3,&[(0,0,1.0),(0,2,2.0),(1,1,3.0)]);
        let mt=m.transpose();assert_eq!(mt.rows,3);assert_eq!(mt.cols,2);assert!(approx_eq(mt.get(2,0),2.0,1e-10));}
    #[test]fn test_cg(){let m=SparseMatrix::from_triplets(2,2,&[(0,0,4.0),(0,1,1.0),(1,0,1.0),(1,1,3.0)]);
        let x=m.solve_cg(&[1.0,2.0],1e-10,100);let ax=m.mul_vec(&x);assert!(approx_eq(ax[0],1.0,1e-6));assert!(approx_eq(ax[1],2.0,1e-6));}
    #[test]fn test_identity(){let id=sparse_identity(3);assert!(approx_eq(id.get(1,1),1.0,1e-10));assert!(approx_eq(id.get(0,1),0.0,1e-10));}
}
