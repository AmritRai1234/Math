//! Set theory: set operations, power sets, relations, cardinality.

use std::collections::HashSet;

/// Power set of a slice (all subsets). Returns 2^n subsets.
pub fn power_set<T:Clone+Eq+std::hash::Hash>(items:&[T])->Vec<Vec<T>>{
    let n=items.len();let mut result=Vec::with_capacity(1<<n);
    for mask in 0..(1u64<<n){let subset:Vec<T>=(0..n).filter(|&i|mask&(1<<i)!=0).map(|i|items[i].clone()).collect();
        result.push(subset);}result
}

/// Union of two sets.
pub fn union<T:Clone+Eq+std::hash::Hash>(a:&[T],b:&[T])->Vec<T>{
    let mut set:HashSet<T>=a.iter().cloned().collect();for x in b{set.insert(x.clone());}set.into_iter().collect()
}

/// Intersection.
pub fn intersection<T:Clone+Eq+std::hash::Hash>(a:&[T],b:&[T])->Vec<T>{
    let sb:HashSet<T>=b.iter().cloned().collect();a.iter().filter(|x|sb.contains(x)).cloned().collect()
}

/// Difference A \ B.
pub fn difference<T:Clone+Eq+std::hash::Hash>(a:&[T],b:&[T])->Vec<T>{
    let sb:HashSet<T>=b.iter().cloned().collect();a.iter().filter(|x|!sb.contains(x)).cloned().collect()
}

/// Symmetric difference.
pub fn symmetric_difference<T:Clone+Eq+std::hash::Hash>(a:&[T],b:&[T])->Vec<T>{
    let mut result=difference(a,b);result.extend(difference(b,a));result
}

/// Cartesian product A × B.
pub fn cartesian_product<T:Clone,U:Clone>(a:&[T],b:&[U])->Vec<(T,U)>{
    a.iter().flat_map(|x|b.iter().map(move|y|(x.clone(),y.clone()))).collect()
}

/// Is subset: A ⊆ B?
pub fn is_subset<T:Eq+std::hash::Hash>(a:&[T],b:&[T])->bool{
    let sb:HashSet<&T>=b.iter().collect();a.iter().all(|x|sb.contains(x))
}

/// Is proper subset: A ⊂ B?
pub fn is_proper_subset<T:Eq+std::hash::Hash>(a:&[T],b:&[T])->bool{a.len()<b.len()&&is_subset(a,b)}

/// Are disjoint?
pub fn are_disjoint<T:Eq+std::hash::Hash>(a:&[T],b:&[T])->bool{
    let sb:HashSet<&T>=b.iter().collect();a.iter().all(|x|!sb.contains(x))
}

/// Jaccard index: |A∩B| / |A∪B|.
pub fn jaccard_index<T:Clone+Eq+std::hash::Hash>(a:&[T],b:&[T])->f64{
    let inter=intersection(a,b).len()as f64;let uni=union(a,b).len()as f64;
    if uni==0.0{1.0}else{inter/uni}
}

/// Partition check: are the sets a valid partition of universal set?
pub fn is_partition<T:Clone+Eq+std::hash::Hash>(sets:&[Vec<T>],universal:&[T])->bool{
    // Union must equal universal, pairwise disjoint
    let mut all:HashSet<T>=HashSet::new();
    for(i,s) in sets.iter().enumerate(){for x in s{if !all.insert(x.clone()){return false;}}// duplicate
        for j in i+1..sets.len(){if !are_disjoint(s,&sets[j]){return false;}}}
    let uni:HashSet<T>=universal.iter().cloned().collect();all==uni
}

/// Equivalence classes from relation (pairs of equivalent elements).
pub fn equivalence_classes<T:Clone+Eq+std::hash::Hash>(elements:&[T],equivalent:&[(T,T)])->Vec<Vec<T>>{
    let mut parent:std::collections::HashMap<usize,usize>=(0..elements.len()).map(|i|(i,i)).collect();
    fn find(parent:&mut std::collections::HashMap<usize,usize>,x:usize)->usize{
        let p=parent[&x];if p!=x{let r=find(parent,p);parent.insert(x,r);r}else{x}}
    let idx=|e:&T|elements.iter().position(|x|x==e);
    for(a,b) in equivalent{if let(Some(ia),Some(ib))=(idx(a),idx(b)){
        let(ra,rb)=(find(&mut parent,ia),find(&mut parent,ib));if ra!=rb{parent.insert(ra,rb);}}}
    let mut classes:std::collections::HashMap<usize,Vec<T>>=std::collections::HashMap::new();
    for i in 0..elements.len(){let r=find(&mut parent,i);classes.entry(r).or_default().push(elements[i].clone());}
    classes.into_values().collect()
}

/// Bell number B(n): number of partitions of an n-element set.
pub fn bell_number(n:usize)->u64{
    if n==0{return 1;}let mut tri=vec![vec![0u64;n+1];n+1];tri[0][0]=1;
    for i in 1..=n{tri[i][0]=tri[i-1][i-1];for j in 1..=i{tri[i][j]=tri[i][j-1]+tri[i-1][j-1];}}tri[n][0]
}

/// Stirling number of the second kind S(n,k): ways to partition n elements into k non-empty subsets.
pub fn stirling_second(n:u64,k:u64)->u64{
    if k==0{return if n==0{1}else{0};}if k>n{return 0;}if k==1||k==n{return 1;}
    k*stirling_second(n-1,k)+stirling_second(n-1,k-1)
}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_power_set(){assert_eq!(power_set(&[1,2,3]).len(),8);}
    #[test]fn test_intersection(){let r=intersection(&[1,2,3],&[2,3,4]);assert_eq!(r,vec![2,3]);}
    #[test]fn test_subset(){assert!(is_subset(&[1,2],&[1,2,3]));assert!(!is_subset(&[1,4],&[1,2,3]));}
    #[test]fn test_bell(){assert_eq!(bell_number(4),15);}
    #[test]fn test_stirling(){assert_eq!(stirling_second(4,2),7);}
    #[test]fn test_jaccard(){assert!((jaccard_index(&[1,2,3],&[2,3,4])-0.5).abs()<1e-10);}
}
