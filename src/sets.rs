//! Set theory operations.

use std::collections::HashSet;
use std::hash::Hash;

/// Set union.
pub fn union<T:Clone+Eq+Hash>(a:&[T],b:&[T])->Vec<T>{
    let mut set:HashSet<_>=a.iter().cloned().collect();for x in b{set.insert(x.clone());}set.into_iter().collect()
}

/// Set intersection.
pub fn intersection<T:Clone+Eq+Hash>(a:&[T],b:&[T])->Vec<T>{
    let sb:HashSet<_>=b.iter().cloned().collect();a.iter().filter(|x|sb.contains(x)).cloned().collect()
}

/// Set difference A \ B.
pub fn difference<T:Clone+Eq+Hash>(a:&[T],b:&[T])->Vec<T>{
    let sb:HashSet<_>=b.iter().cloned().collect();a.iter().filter(|x|!sb.contains(x)).cloned().collect()
}

/// Symmetric difference.
pub fn symmetric_diff<T:Clone+Eq+Hash>(a:&[T],b:&[T])->Vec<T>{
    let sa:HashSet<_>=a.iter().cloned().collect();let sb:HashSet<_>=b.iter().cloned().collect();
    sa.symmetric_difference(&sb).cloned().collect()
}

/// Subset check.
pub fn is_subset<T:Clone+Eq+Hash>(a:&[T],b:&[T])->bool{
    let sb:HashSet<_>=b.iter().cloned().collect();a.iter().all(|x|sb.contains(x))
}

/// Power set.
pub fn power_set<T:Clone>(s:&[T])->Vec<Vec<T>>{
    let n=s.len();let count=1<<n;let mut result=Vec::with_capacity(count);
    for mask in 0..count{let mut subset=Vec::new();for i in 0..n{if mask&(1<<i)!=0{subset.push(s[i].clone());}}result.push(subset);}
    result
}

/// Cartesian product.
pub fn cartesian_product<T:Clone,U:Clone>(a:&[T],b:&[U])->Vec<(T,U)>{
    a.iter().flat_map(|x|b.iter().map(move|y|(x.clone(),y.clone()))).collect()
}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_union(){let mut r=union(&[1,2,3],&[3,4,5]);r.sort();assert_eq!(r,vec![1,2,3,4,5]);}
    #[test]fn test_intersection(){let mut r=intersection(&[1,2,3],&[2,3,4]);r.sort();assert_eq!(r,vec![2,3]);}
    #[test]fn test_difference(){let r=difference(&[1,2,3],&[2,3,4]);assert_eq!(r,vec![1]);}
    #[test]fn test_subset(){assert!(is_subset(&[1,2],&[1,2,3]));assert!(!is_subset(&[1,4],&[1,2,3]));}
    #[test]fn test_power_set(){assert_eq!(power_set(&[1,2]).len(),4);}
}
