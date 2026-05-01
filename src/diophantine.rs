//! Diophantine equations: Pythagorean triples, sum of squares, Pell, Chinese Remainder.

/// Generate all primitive Pythagorean triples with hypotenuse ≤ limit.
pub fn pythagorean_triples(limit:u64)->Vec<(u64,u64,u64)>{
    let mut triples=Vec::new();let mut m=2u64;
    while m*m<limit{for n in 1..m{if(m-n)%2==1&&gcd(m,n)==1{
        let(a,b,c)=(m*m-n*n,2*m*n,m*m+n*n);if c<=limit{triples.push(if a<b{(a,b,c)}else{(b,a,c)});}}}m+=1;}
    triples.sort();triples
}

fn gcd(mut a:u64,mut b:u64)->u64{while b!=0{let t=b;b=a%b;a=t;}a}

/// Check if n is sum of two squares.
pub fn is_sum_of_two_squares(n:u64)->bool{
    let mut i=0u64;while i*i<=n{let rem=n-i*i;
        let s=crate::powers::sqrt(rem as f64)as u64;if s*s==rem{return true;}i+=1;}false
}

/// Find representation n = a² + b² if it exists.
pub fn sum_of_two_squares(n:u64)->Option<(u64,u64)>{
    let mut i=0u64;while i*i<=n{let rem=n-i*i;
        let s=crate::powers::sqrt(rem as f64)as u64;if s*s==rem{return Some((i,s));}i+=1;}None
}

/// Lagrange's four-square theorem: represent n as sum of 4 squares.
pub fn four_squares(n:u64)->(u64,u64,u64,u64){
    let isq=crate::powers::sqrt(n as f64)as u64;
    for a in (0..=isq).rev(){let r1=n-a*a;
        let bmax=crate::powers::sqrt(r1 as f64)as u64;
        for b in (0..=bmax).rev(){let r2=r1-b*b;
            if let Some((c,d))=sum_of_two_squares(r2){return(a,b,c,d);}}}
    (0,0,0,0)
}

/// Chinese Remainder Theorem: solve x ≡ a_i (mod m_i) for pairwise coprime m_i.
pub fn chinese_remainder(residues:&[i64],moduli:&[i64])->Option<i64>{
    let n=residues.len();if n==0{return None;}
    let prod:i64=moduli.iter().product();
    let mut sum=0i64;
    for i in 0..n{let mi=moduli[i];let ni=prod/mi;
        let(_,xi,_)=extended_gcd(ni,mi);sum=(sum+residues[i]*ni%prod*xi%prod)%prod;}
    Some(((sum%prod)+prod)%prod)
}

fn extended_gcd(a:i64,b:i64)->(i64,i64,i64){
    if a==0{return(b,0,1);}let(g,x,y)=extended_gcd(b%a,a);(g,y-(b/a)*x,x)
}

/// Solve linear Diophantine equation ax + by = c.
pub fn linear_diophantine(a:i64,b:i64,c:i64)->Option<(i64,i64)>{
    let(g,x0,y0)=extended_gcd(a.abs(),b.abs());
    if c%g!=0{return None;}let scale=c/g;
    let x=x0*scale*(if a<0{-1}else{1});let y=y0*scale*(if b<0{-1}else{1});Some((x,y))
}

/// Check if number is a perfect power (n = a^k for some k ≥ 2).
pub fn is_perfect_power(n:u64)->bool{
    for k in 2..=63{let root=crate::powers::pow(n as f64,1.0/k as f64).round()as u64;
        if crate::number_theory::mod_pow(root,k,u64::MAX)==n{return true;}
        if root>1&&{let mut p=1u64;let mut ok=true;for _ in 0..k{
            p=match p.checked_mul(root){Some(v)=>v,None=>{ok=false;break;0}};
            if !ok{break;}}ok&&p==n}{return true;}}false
}

/// Fermat factorization.
pub fn fermat_factor(n:u64)->Option<(u64,u64)>{
    if n%2==0{return Some((2,n/2));}
    let mut a=crate::powers::sqrt(n as f64).ceil()as u64;
    for _ in 0..10000{let b2=a*a-n;let b=crate::powers::sqrt(b2 as f64)as u64;
        if b*b==b2{return Some((a-b,a+b));}a+=1;}None
}

/// Count integer points in circle of radius r: |{(x,y) : x² + y² ≤ r²}|.
pub fn lattice_points_in_circle(r:u64)->u64{
    let mut count=0u64;let r2=r*r;
    for x in 0..=r{let y_max=crate::powers::sqrt((r2-x*x)as f64)as u64;
        count+=y_max+1;}// first quadrant + axes
    4*(count-r-1)+4*r+1 // all quadrants
}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_pythagorean(){let t=pythagorean_triples(100);assert!(t.contains(&(3,4,5)));}
    #[test]fn test_sum_squares(){assert!(is_sum_of_two_squares(25));assert_eq!(sum_of_two_squares(25),Some((0,5)));}
    #[test]fn test_four_squares(){let(a,b,c,d)=four_squares(15);assert_eq!(a*a+b*b+c*c+d*d,15);}
    #[test]fn test_crt(){assert_eq!(chinese_remainder(&[2,3,2],&[3,5,7]),Some(23));}
    #[test]fn test_diophantine(){let(x,y)=linear_diophantine(3,5,1).unwrap();assert_eq!(3*x+5*y,1);}
}
