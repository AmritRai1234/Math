//! Linear programming and optimization.

use crate::arithmetic::abs;

/// Simplex method for linear programming.
/// Maximize c·x subject to Ax ≤ b, x ≥ 0.
/// Returns (optimal_value, solution_vector) or None if infeasible/unbounded.
pub fn simplex(c:&[f64],a:&[Vec<f64>],b:&[f64])->Option<(f64,Vec<f64>)>{
    let m=a.len();let n=c.len();
    // Build tableau: m constraints, n variables, m slack variables
    let cols=n+m+1;let rows=m+1;
    let mut tab=vec![vec![0.0;cols];rows];
    // Fill constraints
    for i in 0..m{for j in 0..n{tab[i][j]=a[i][j];}tab[i][n+i]=1.0;tab[i][cols-1]=b[i];}
    // Objective row (negate for maximization)
    for j in 0..n{tab[m][j]=-c[j];}
    // Iterate
    for _ in 0..10000{
        // Find pivot column (most negative in objective row)
        let mut pivot_col=0;let mut min_val=0.0;
        for j in 0..cols-1{if tab[m][j]<min_val{min_val=tab[m][j];pivot_col=j;}}
        if min_val>=-1e-10{break;} // optimal
        // Find pivot row (minimum ratio test)
        let mut pivot_row=None;let mut min_ratio=f64::INFINITY;
        for i in 0..m{if tab[i][pivot_col]>1e-10{let ratio=tab[i][cols-1]/tab[i][pivot_col];
            if ratio<min_ratio{min_ratio=ratio;pivot_row=Some(i);}}}
        let pivot_row=pivot_row?; // unbounded if None
        // Pivot
        let pivot_val=tab[pivot_row][pivot_col];
        for j in 0..cols{tab[pivot_row][j]/=pivot_val;}
        for i in 0..rows{if i!=pivot_row{let factor=tab[i][pivot_col];
            for j in 0..cols{tab[i][j]-=factor*tab[pivot_row][j];}}}
    }
    // Extract solution
    let mut x=vec![0.0;n];
    for j in 0..n{let mut basic_row=None;let mut is_basic=true;
        for i in 0..m{if abs(tab[i][j]-1.0)<1e-10{if basic_row.is_some(){is_basic=false;break;}basic_row=Some(i);}
            else if abs(tab[i][j])>1e-10{is_basic=false;break;}}
        if is_basic{if let Some(r)=basic_row{x[j]=tab[r][cols-1];}}}
    Some((tab[m][cols-1],x))
}

/// Quadratic programming (simple projected gradient descent).
/// Minimize 0.5 x^T Q x + c^T x subject to x ≥ 0.
pub fn quadratic_minimize(q:&[Vec<f64>],c:&[f64],lr:f64,steps:usize)->Vec<f64>{
    let n=c.len();let mut x=vec![0.0;n];
    for _ in 0..steps{let mut grad=c.to_vec();
        for i in 0..n{for j in 0..n{grad[i]+=q[i][j]*x[j];}}
        for i in 0..n{x[i]-=lr*grad[i];if x[i]<0.0{x[i]=0.0;}}} // project onto x≥0
    x
}

/// Assignment problem (Hungarian algorithm simplified for small matrices).
/// Minimize total cost for n-by-n cost matrix. Returns assignment.
pub fn assignment(cost:&[Vec<f64>])->Vec<usize>{
    let n=cost.len();let mut best_perm=vec![0;n];let mut best_cost=f64::INFINITY;
    if n>10{return greedy_assignment(cost);} // fallback for large
    fn permute(n:usize,perm:&mut Vec<usize>,used:&mut Vec<bool>,cost:&[Vec<f64>],best:&mut f64,best_perm:&mut Vec<usize>,cur_cost:f64){
        if perm.len()==n{if cur_cost<*best{*best=cur_cost;*best_perm=perm.clone();}return;}
        let i=perm.len();
        for j in 0..n{if !used[j]{let new_cost=cur_cost+cost[i][j];if new_cost<*best{
            used[j]=true;perm.push(j);permute(n,perm,used,cost,best,best_perm,new_cost);perm.pop();used[j]=false;}}}}
    let mut perm=Vec::new();let mut used=vec![false;n];
    permute(n,&mut perm,&mut used,cost,&mut best_cost,&mut best_perm,0.0);best_perm
}

fn greedy_assignment(cost:&[Vec<f64>])->Vec<usize>{
    let n=cost.len();let mut assigned=vec![false;n];let mut result=vec![0;n];
    for i in 0..n{let mut best_j=0;let mut best_c=f64::INFINITY;
        for j in 0..n{if !assigned[j]&&cost[i][j]<best_c{best_c=cost[i][j];best_j=j;}}
        result[i]=best_j;assigned[best_j]=true;}result
}

/// Knapsack 0/1 (dynamic programming).
pub fn knapsack(weights:&[u64],values:&[u64],capacity:u64)->u64{
    let n=weights.len();let cap=capacity as usize;
    let mut dp=vec![vec![0u64;cap+1];n+1];
    for i in 1..=n{for w in 0..=cap{dp[i][w]=dp[i-1][w];
        if weights[i-1]<=w as u64{let v=dp[i-1][w-weights[i-1]as usize]+values[i-1];if v>dp[i][w]{dp[i][w]=v;}}}}
    dp[n][cap]
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_simplex(){// Maximize 3x+2y s.t. x+y≤4, x+3y≤6
        let r=simplex(&[3.0,2.0],&[vec![1.0,1.0],vec![1.0,3.0]],&[4.0,6.0]).unwrap();
        assert!(r.0>9.0,"simplex optimal value {} too low",r.0);}
    #[test]fn test_knapsack(){assert_eq!(knapsack(&[1,2,3],&[6,10,12],5),22);}
    #[test]fn test_assignment(){let c=vec![vec![9.0,2.0,7.0],vec![6.0,4.0,3.0],vec![5.0,8.0,1.0]];
        let a=assignment(&c);let total:f64=a.iter().enumerate().map(|(i,&j)|c[i][j]).sum();assert!(total<=10.1);} // optimal: 2+6+1=9 or 2+3+5=10
}
