//! Abstract algebra: groups, rings, permutations.

/// Permutation represented as a mapping.
#[derive(Clone,Debug,PartialEq)]
pub struct Permutation{pub map:Vec<usize>}

impl Permutation{
    pub fn identity(n:usize)->Self{Self{map:(0..n).collect()}}
    pub fn from_vec(v:Vec<usize>)->Self{Self{map:v}}
    pub fn len(&self)->usize{self.map.len()}
    pub fn apply(&self,i:usize)->usize{self.map[i]}

    /// Compose: self ∘ other (apply other first, then self).
    pub fn compose(&self,other:&Self)->Self{Self{map:other.map.iter().map(|&i|self.map[i]).collect()}}

    /// Inverse permutation.
    pub fn inverse(&self)->Self{let mut inv=vec![0;self.len()];for(i,&j) in self.map.iter().enumerate(){inv[j]=i;}Self{map:inv}}

    /// Cycle decomposition.
    pub fn cycles(&self)->Vec<Vec<usize>>{
        let n=self.len();let mut visited=vec![false;n];let mut cycles=Vec::new();
        for i in 0..n{if visited[i]{continue;}let mut cycle=vec![i];visited[i]=true;
            let mut j=self.map[i];while j!=i{cycle.push(j);visited[j]=true;j=self.map[j];}
            if cycle.len()>1{cycles.push(cycle);}}cycles
    }

    /// Order (smallest k such that σ^k = identity).
    pub fn order(&self)->usize{
        let cycles=self.cycles();if cycles.is_empty(){return 1;}
        cycles.iter().map(|c|c.len()).fold(1,lcm)
    }

    /// Sign: +1 for even, -1 for odd.
    pub fn sign(&self)->i32{let n_transpositions:usize=self.cycles().iter().map(|c|c.len()-1).sum();if n_transpositions%2==0{1}else{-1}}

    /// Check if this is a valid permutation.
    pub fn is_valid(&self)->bool{let n=self.len();let mut seen=vec![false;n];
        for &x in &self.map{if x>=n||seen[x]{return false;}seen[x]=true;}true}

    /// Apply permutation to a slice.
    pub fn permute<T:Clone>(&self,data:&[T])->Vec<T>{self.map.iter().map(|&i|data[i].clone()).collect()}

    /// Transposition (swap i and j).
    pub fn transposition(n:usize,i:usize,j:usize)->Self{let mut map:Vec<usize>=(0..n).collect();map[i]=j;map[j]=i;Self{map}}

    /// Power: σ^k.
    pub fn pow(&self,k:u32)->Self{let mut result=Self::identity(self.len());let mut base=self.clone();
        let mut exp=k;while exp>0{if exp&1==1{result=result.compose(&base);}base=base.compose(&base);exp>>=1;}result}

    /// Number of inversions.
    pub fn inversions(&self)->usize{let n=self.len();let mut count=0;
        for i in 0..n{for j in i+1..n{if self.map[i]>self.map[j]{count+=1;}}}count}
}

fn gcd(a:usize,b:usize)->usize{if b==0{a}else{gcd(b,a%b)}}
fn lcm(a:usize,b:usize)->usize{a/gcd(a,b)*b}

/// Generate all permutations of n elements.
pub fn all_permutations(n:usize)->Vec<Permutation>{
    let mut result=Vec::new();let mut perm:Vec<usize>=(0..n).collect();
    fn gen(perm:&mut Vec<usize>,k:usize,result:&mut Vec<Permutation>){
        if k==1{result.push(Permutation::from_vec(perm.clone()));return;}
        for i in 0..k{gen(perm,k-1,result);if k%2==0{perm.swap(i,k-1);}else{perm.swap(0,k-1);}}}
    gen(&mut perm,n,&mut result);result
}

/// Symmetric group S_n order = n!
pub fn symmetric_group_order(n:u64)->u64{crate::combinatorics::factorial(n)}

/// Check if a set with an operation forms a group (closure, assoc, identity, inverse).
/// Elements are indices 0..n, cayley_table[i][j] = i*j.
pub fn is_group(cayley_table:&[Vec<usize>])->bool{
    let n=cayley_table.len();
    // Closure (all entries < n)
    for row in cayley_table{for &x in row{if x>=n{return false;}}}
    // Identity exists
    let mut identity=None;
    for e in 0..n{if(0..n).all(|i|cayley_table[e][i]==i&&cayley_table[i][e]==i){identity=Some(e);break;}}
    let e=match identity{Some(e)=>e,None=>return false};
    // Inverse for each element
    for i in 0..n{if !(0..n).any(|j|cayley_table[i][j]==e&&cayley_table[j][i]==e){return false;}}
    true
}

/// Check if group is abelian (commutative).
pub fn is_abelian(cayley_table:&[Vec<usize>])->bool{
    let n=cayley_table.len();
    for i in 0..n{for j in 0..n{if cayley_table[i][j]!=cayley_table[j][i]{return false;}}}true
}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_identity(){let p=Permutation::identity(4);assert_eq!(p.apply(2),2);}
    #[test]fn test_compose(){let p=Permutation::from_vec(vec![1,0,2]);let q=Permutation::from_vec(vec![0,2,1]);
        let r=p.compose(&q);assert_eq!(r.map,vec![1,2,0]);}
    #[test]fn test_inverse(){let p=Permutation::from_vec(vec![1,2,0]);let inv=p.inverse();
        let id=p.compose(&inv);assert_eq!(id,Permutation::identity(3));}
    #[test]fn test_cycles(){let p=Permutation::from_vec(vec![1,2,0,3]);let c=p.cycles();assert_eq!(c.len(),1);assert_eq!(c[0].len(),3);}
    #[test]fn test_sign(){let p=Permutation::transposition(3,0,1);assert_eq!(p.sign(),-1);}
    #[test]fn test_order(){let p=Permutation::from_vec(vec![1,2,0]);assert_eq!(p.order(),3);}
    #[test]fn test_is_group(){// Z/2Z: {0,1} with XOR
        let table=vec![vec![0,1],vec![1,0]];assert!(is_group(&table));assert!(is_abelian(&table));}
}
