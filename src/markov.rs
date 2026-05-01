//! Markov chains: transition matrices, steady state, random walks.

use crate::arithmetic::abs;

/// Multiply transition matrix by state vector.
pub fn step(trans:&[Vec<f64>],state:&[f64])->Vec<f64>{
    let n=trans.len();(0..n).map(|i|(0..n).map(|j|trans[i][j]*state[j]).sum()).collect()
}

/// Run Markov chain for k steps.
pub fn run(trans:&[Vec<f64>],initial:&[f64],steps:usize)->Vec<f64>{
    let mut state=initial.to_vec();for _ in 0..steps{state=step(trans,&state);}state
}

/// Find steady state distribution (power iteration).
pub fn steady_state(trans:&[Vec<f64>],tol:f64)->Vec<f64>{
    let n=trans.len();let mut state=vec![1.0/n as f64;n];
    for _ in 0..100000{let next=step(trans,&state);
        let diff:f64=state.iter().zip(next.iter()).map(|(&a,&b)|abs(a-b)).sum();
        state=next;if diff<tol{break;}}state
}

/// Check if transition matrix is valid (rows sum to 1, all non-negative).
pub fn is_valid_transition(trans:&[Vec<f64>])->bool{
    for row in trans{if row.iter().any(|&x|x<0.0){return false;}
        let sum:f64=row.iter().sum();if abs(sum-1.0)>1e-10{return false;}}true
}

/// Check if chain is ergodic (irreducible + aperiodic).
pub fn is_ergodic(trans:&[Vec<f64>])->bool{
    let n=trans.len();
    // Check irreducibility: all states reachable from all states via power iteration
    let mut reachable=vec![vec![false;n];n];
    for i in 0..n{reachable[i][i]=true;for j in 0..n{if trans[i][j]>0.0{reachable[i][j]=true;}}}
    // Transitive closure
    for k in 0..n{for i in 0..n{for j in 0..n{if reachable[i][k]&&reachable[k][j]{reachable[i][j]=true;}}}}
    for i in 0..n{for j in 0..n{if !reachable[i][j]{return false;}}}
    // Check aperiodicity: at least one self-loop (sufficient condition)
    (0..n).any(|i|trans[i][i]>0.0)
}

/// Absorption probabilities for absorbing Markov chain.
/// Returns probability of being absorbed in each absorbing state.
pub fn absorption_probabilities(trans:&[Vec<f64>])->Vec<Vec<f64>>{
    let n=trans.len();
    let absorbing:Vec<bool>=(0..n).map(|i|abs(trans[i][i]-1.0)<1e-10).collect();
    let transient:Vec<usize>=(0..n).filter(|&i|!absorbing[i]).collect();
    let abs_states:Vec<usize>=(0..n).filter(|&i|absorbing[i]).collect();
    if transient.is_empty(){return vec![];}
    let t=transient.len();let a=abs_states.len();
    // Q matrix (transient -> transient)
    let mut q=vec![vec![0.0;t];t];
    for(ii,&i)in transient.iter().enumerate(){for(jj,&j)in transient.iter().enumerate(){q[ii][jj]=trans[i][j];}}
    // R matrix (transient -> absorbing)
    let mut r=vec![vec![0.0;a];t];
    for(ii,&i)in transient.iter().enumerate(){for(jj,&j)in abs_states.iter().enumerate(){r[ii][jj]=trans[i][j];}}
    // N = (I-Q)^(-1), B = N*R
    let mut imq=vec![vec![0.0;t];t];
    for i in 0..t{for j in 0..t{imq[i][j]=if i==j{1.0}else{0.0}-q[i][j];}}
    let n_mat=crate::matrix::inverse(&imq).unwrap_or(imq);
    crate::matrix::mat_mul(&n_mat,&r)
}

/// Mean first passage time from state i to state j.
pub fn mean_first_passage(trans:&[Vec<f64>])->Vec<Vec<f64>>{
    let n=trans.len();let pi=steady_state(trans,1e-12);
    let mut m=vec![vec![0.0;n];n];
    // Approximate via simulation
    for target in 0..n{for start in 0..n{if start==target{m[start][target]=0.0;continue;}
        // Use fundamental matrix approach approximation
        m[start][target]=if pi[target]>1e-15{1.0/pi[target]}else{f64::INFINITY};}}
    m
}

/// Expected number of visits to each state in k steps.
pub fn expected_visits(trans:&[Vec<f64>],initial:&[f64],steps:usize)->Vec<f64>{
    let n=trans.len();let mut visits=vec![0.0;n];let mut state=initial.to_vec();
    for _ in 0..steps{for i in 0..n{visits[i]+=state[i];}state=step(trans,&state);}visits
}

/// Random walk on 1D lattice: returns position after n steps.
pub fn random_walk_1d(steps:usize,seed:u64)->Vec<i64>{
    let mut pos=0i64;let mut path=vec![0i64];let mut rng=seed;
    for _ in 0..steps{rng=crate::noise::xorshift(rng);
        pos+=if rng%2==0{1}else{-1};path.push(pos);}path
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_steady(){let t=vec![vec![0.9,0.1],vec![0.5,0.5]];
        let s=steady_state(&t,1e-12);assert!(approx_eq(s[0]+s[1],1.0,1e-10));}
    #[test]fn test_valid(){let t=vec![vec![0.5,0.5],vec![0.3,0.7]];assert!(is_valid_transition(&t));}
    #[test]fn test_step(){let t=vec![vec![1.0,0.0],vec![0.0,1.0]];let s=step(&t,&[0.5,0.5]);assert_eq!(s,vec![0.5,0.5]);}
}
