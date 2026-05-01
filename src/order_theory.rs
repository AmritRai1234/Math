//! Order theory: lattices, partial orders, Möbius function.

use std::collections::HashSet;

/// Partial order represented as a relation matrix (relation[i][j] = true means i ≤ j).
pub type Relation=Vec<Vec<bool>>;

/// Check if relation is a partial order (reflexive, antisymmetric, transitive).
pub fn is_partial_order(r:&Relation)->bool{
    let n=r.len();
    // Reflexive
    for i in 0..n{if !r[i][i]{return false;}}
    // Antisymmetric
    for i in 0..n{for j in 0..n{if i!=j&&r[i][j]&&r[j][i]{return false;}}}
    // Transitive
    for i in 0..n{for j in 0..n{for k in 0..n{if r[i][j]&&r[j][k]&&!r[i][k]{return false;}}}}
    true
}

/// Check if relation is a total order.
pub fn is_total_order(r:&Relation)->bool{
    if !is_partial_order(r){return false;}
    let n=r.len();for i in 0..n{for j in 0..n{if !r[i][j]&&!r[j][i]{return false;}}}true
}

/// Transitive closure of a relation.
pub fn transitive_closure(r:&Relation)->Relation{
    let n=r.len();let mut result=r.clone();
    for k in 0..n{for i in 0..n{for j in 0..n{if result[i][k]&&result[k][j]{result[i][j]=true;}}}}result
}

/// Hasse diagram: remove transitive edges from partial order.
pub fn hasse_diagram(r:&Relation)->Relation{
    let n=r.len();let mut h=r.clone();
    for i in 0..n{h[i][i]=false;} // remove reflexive
    for i in 0..n{for j in 0..n{if h[i][j]{for k in 0..n{if h[i][k]&&h[k][j]&&k!=i&&k!=j{h[i][j]=false;}}}}}
    h
}

/// Join (least upper bound) of two elements. Returns None if doesn't exist.
pub fn join(r:&Relation,a:usize,b:usize)->Option<usize>{
    let n=r.len();
    let upper_bounds:Vec<usize>=(0..n).filter(|&c|r[a][c]&&r[b][c]).collect();
    upper_bounds.iter().find(|&&c|upper_bounds.iter().all(|&d|r[c][d])).cloned()
}

/// Meet (greatest lower bound).
pub fn meet(r:&Relation,a:usize,b:usize)->Option<usize>{
    let n=r.len();
    let lower_bounds:Vec<usize>=(0..n).filter(|&c|r[c][a]&&r[c][b]).collect();
    lower_bounds.iter().find(|&&c|lower_bounds.iter().all(|&d|r[d][c])).cloned()
}

/// Check if poset is a lattice (every pair has join and meet).
pub fn is_lattice(r:&Relation)->bool{
    let n=r.len();
    for i in 0..n{for j in 0..n{if join(r,i,j).is_none()||meet(r,i,j).is_none(){return false;}}}true
}

/// Minimal elements.
pub fn minimal_elements(r:&Relation)->Vec<usize>{
    let n=r.len();(0..n).filter(|&i|(0..n).all(|j|!r[j][i]||j==i)).collect()
}

/// Maximal elements.
pub fn maximal_elements(r:&Relation)->Vec<usize>{
    let n=r.len();(0..n).filter(|&i|(0..n).all(|j|!r[i][j]||j==i)).collect()
}

/// Möbius function μ(x,y) for a poset.
pub fn mobius_function(r:&Relation)->Vec<Vec<i64>>{
    let n=r.len();let mut mu=vec![vec![0i64;n];n];
    for i in 0..n{mu[i][i]=1;
        for j in 0..n{if i!=j&&r[i][j]{let between:Vec<usize>=(0..n).filter(|&k|k!=j&&r[i][k]&&r[k][j]).collect();
            mu[i][j]=-between.iter().map(|&k|mu[i][k]).sum::<i64>();}}}mu
}

/// Zeta function ζ(x,y) = 1 if x ≤ y, else 0.
pub fn zeta_function(r:&Relation)->Vec<Vec<i64>>{
    let n=r.len();(0..n).map(|i|(0..n).map(|j|if r[i][j]{1}else{0}).collect()).collect()
}

/// Topological sort of poset.
pub fn topological_sort(r:&Relation)->Vec<usize>{
    let n=r.len();let mut in_degree=vec![0usize;n];
    for j in 0..n{for i in 0..n{if i!=j&&r[i][j]{in_degree[j]+=1;}}}
    let mut result=Vec::new();let mut available:Vec<usize>=(0..n).filter(|&i|in_degree[i]==0).collect();
    while let Some(node)=available.pop(){result.push(node);
        for j in 0..n{if node!=j&&r[node][j]{in_degree[j]-=1;if in_degree[j]==0{available.push(j);}}}}result
}

/// Width of poset (size of largest antichain) — Dilworth's theorem gives min chain decomposition.
pub fn poset_width(r:&Relation)->usize{
    let n=r.len();let mut max_antichain=0;
    // Brute-force for small n
    for mask in 1..(1u64<<n){let elements:Vec<usize>=(0..n).filter(|&i|mask&(1<<i)!=0).collect();
        let is_antichain=elements.iter().all(|&i|elements.iter().all(|&j|i==j||(!r[i][j]&&!r[j][i])));
        if is_antichain{max_antichain=max_antichain.max(elements.len());}}
    max_antichain
}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_divisibility(){// Divisibility on {1,2,3,6}: 1|2, 1|3, 1|6, 2|6, 3|6
        let r=vec![vec![true,true,true,true],vec![false,true,false,true],vec![false,false,true,true],vec![false,false,false,true]];
        assert!(is_partial_order(&r));assert!(is_lattice(&r));}
    #[test]fn test_minimal(){let r=vec![vec![true,true,true],vec![false,true,true],vec![false,false,true]];
        assert_eq!(minimal_elements(&r),vec![0]);}
    #[test]fn test_mobius(){let r=vec![vec![true,true],vec![false,true]];let mu=mobius_function(&r);assert_eq!(mu[0][1],-1);}
}
