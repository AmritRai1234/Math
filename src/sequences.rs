//! Number sequences beyond Fibonacci.

/// Tribonacci: T(0)=0, T(1)=0, T(2)=1, T(n)=T(n-1)+T(n-2)+T(n-3).
pub fn tribonacci(n:u64)->u64{
    if n==0||n==1{return 0;}if n==2{return 1;}
    let(mut a,mut b,mut c)=(0u64,0u64,1u64);
    for _ in 3..=n{let t=a.saturating_add(b).saturating_add(c);a=b;b=c;c=t;}c
}

/// Triangular numbers: T(n) = n*(n+1)/2.
pub fn triangular(n:u64)->u64{n*(n+1)/2}

/// Square numbers.
pub fn square_number(n:u64)->u64{n*n}

/// Pentagonal numbers: P(n) = n*(3n-1)/2.
pub fn pentagonal(n:u64)->u64{n*(3*n-1)/2}

/// Hexagonal numbers: H(n) = n*(2n-1).
pub fn hexagonal(n:u64)->u64{n*(2*n-1)}

/// Tetrahedral numbers: T(n) = n*(n+1)*(n+2)/6.
pub fn tetrahedral(n:u64)->u64{n*(n+1)*(n+2)/6}

/// Pyramidal (square) numbers: P(n) = n*(n+1)*(2n+1)/6.
pub fn pyramidal(n:u64)->u64{n*(n+1)*(2*n+1)/6}

/// Pronic numbers (oblong): n*(n+1).
pub fn pronic(n:u64)->u64{n*(n+1)}

/// Look-and-say sequence (returns nth term as string).
pub fn look_and_say(n:u32)->String{
    let mut s=String::from("1");
    for _ in 1..n{let mut next=String::new();let chars:Vec<char>=s.chars().collect();
        let mut i=0;while i<chars.len(){let c=chars[i];let mut count=1;
            while i+count<chars.len()&&chars[i+count]==c{count+=1;}
            next.push_str(&count.to_string());next.push(c);i+=count;}s=next;}
    s
}

/// Collatz sequence length (number of steps to reach 1).
pub fn collatz_length(mut n:u64)->u64{
    let mut count=0;while n!=1&&n!=0{n=if n%2==0{n/2}else{3*n+1};count+=1;}count
}

/// Stern-Brocot sequence.
pub fn stern_brocot(n:usize)->Vec<u64>{
    if n==0{return vec![];}let mut seq=vec![1u64,1];
    while seq.len()<n{let i=seq.len()/2;
        seq.push(seq[i-1]+seq[i]);if seq.len()<n{seq.push(seq[i]);}}
    seq.truncate(n);seq
}

/// Recaman's sequence.
pub fn recaman(n:usize)->Vec<u64>{
    let mut seq=vec![0u64];let mut seen=std::collections::HashSet::new();seen.insert(0u64);
    for i in 1..n{let prev=seq[i-1];let sub=if prev>i as u64{prev-i as u64}else{0};
        let val=if sub>0&&!seen.contains(&sub){sub}else{prev+i as u64};
        seq.push(val);seen.insert(val);}
    seq
}

/// Sylvester's sequence: a(0)=2, a(n)=a(n-1)*(a(n-1)-1)+1.
pub fn sylvester(n:u32)->u64{
    let mut a=2u64;for _ in 1..=n{a=a.saturating_mul(a-1).saturating_add(1);}a
}

/// Padovan sequence: P(0)=P(1)=P(2)=1, P(n)=P(n-2)+P(n-3).
pub fn padovan(n:u64)->u64{
    if n<=2{return 1;}let(mut a,mut b,mut c)=(1u64,1u64,1u64);
    for _ in 3..=n{let t=a.saturating_add(b);a=b;b=c;c=t;}c
}

/// Pell numbers: P(0)=0, P(1)=1, P(n)=2*P(n-1)+P(n-2).
pub fn pell(n:u64)->u64{
    if n==0{return 0;}if n==1{return 1;}
    let(mut a,mut b)=(0u64,1u64);for _ in 2..=n{let t=2u64.saturating_mul(b).saturating_add(a);a=b;b=t;}b
}

/// Jacobsthal numbers: J(0)=0, J(1)=1, J(n)=J(n-1)+2*J(n-2).
pub fn jacobsthal(n:u64)->u64{
    if n==0{return 0;}if n==1{return 1;}
    let(mut a,mut b)=(0u64,1u64);for _ in 2..=n{let t=b.saturating_add(2u64.saturating_mul(a));a=b;b=t;}b
}

/// Motzkin numbers.
pub fn motzkin(n:u64)->u64{
    if n<=1{return 1;}let n=n as usize;
    let mut m=vec![0u64;n+1];m[0]=1;m[1]=1;
    for i in 2..=n{m[i]=m[i-1]+((0..i-1).map(|k|m[k]*m[i-2-k]).sum::<u64>());}m[n]
}

/// Narayana numbers N(n,k).
pub fn narayana(n:u64,k:u64)->u64{
    if k==0||k>n{return 0;}
    crate::combinatorics::combinations(n,k)*crate::combinatorics::combinations(n,k-1)/(n)
}

/// Harmonic number H(n) = 1 + 1/2 + 1/3 + ... + 1/n.
pub fn harmonic(n:u64)->f64{
    let mut sum=0.0;for i in 1..=n{sum+=1.0/i as f64;}sum
}

/// Bernoulli numbers B(n) (first few).
pub fn bernoulli(n:u32)->f64{
    let mut a=vec![0.0f64;(n+1)as usize];
    for m in 0..=n as usize{a[m]=1.0/((m+1)as f64);
        for j in(1..=m).rev(){a[j-1]=(j as f64)*(a[j-1]-a[j]);}}
    a[0]
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_tribonacci(){assert_eq!(tribonacci(0),0);assert_eq!(tribonacci(4),2);assert_eq!(tribonacci(7),13);}
    #[test]fn test_triangular(){assert_eq!(triangular(4),10);assert_eq!(triangular(100),5050);}
    #[test]fn test_pell(){assert_eq!(pell(0),0);assert_eq!(pell(1),1);assert_eq!(pell(5),29);}
    #[test]fn test_padovan(){assert_eq!(padovan(0),1);assert_eq!(padovan(7),5);}
    #[test]fn test_look_and_say(){assert_eq!(look_and_say(1),"1");assert_eq!(look_and_say(2),"11");assert_eq!(look_and_say(3),"21");}
    #[test]fn test_harmonic(){assert!(approx_eq(harmonic(1),1.0,1e-14));assert!(harmonic(10)>2.9&&harmonic(10)<3.0);}
    #[test]fn test_bernoulli(){assert!(approx_eq(bernoulli(0),1.0,1e-10));assert!(approx_eq(bernoulli(1),0.5,1e-10));assert!(approx_eq(bernoulli(2),1.0/6.0,1e-10));}
}
