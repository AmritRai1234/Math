//! Continued fractions: convergents, best rational approximations.

/// Regular continued fraction expansion of a real number.
pub fn continued_fraction(x:f64,max_terms:usize)->Vec<i64>{
    let mut result=Vec::new();let mut x=x;
    for _ in 0..max_terms{let a=x.floor()as i64;result.push(a);
        let frac=x-a as f64;if frac.abs()<1e-12{break;}x=1.0/frac;}result
}

/// Evaluate continued fraction [a0; a1, a2, ...] to f64.
pub fn evaluate_cf(cf:&[i64])->f64{
    let mut result=0.0;for &a in cf.iter().rev(){result=1.0/(result+a as f64);}
    1.0/result
}

/// Convergents (best rational approximations) from continued fraction.
pub fn convergents(cf:&[i64])->Vec<(i64,i64)>{
    let mut result=Vec::new();let(mut h0,mut h1)=(0i64,1i64);let(mut k0,mut k1)=(1i64,0i64);
    for &a in cf{let h=a*h1+h0;let k=a*k1+k0;result.push((h,k));h0=h1;h1=h;k0=k1;k1=k;}result
}

/// Best rational approximation of x with denominator ≤ max_denom.
pub fn best_rational(x:f64,max_denom:i64)->(i64,i64){
    let cf=continued_fraction(x,30);let convs=convergents(&cf);
    let mut best=(x.round()as i64,1i64);
    for &(h,k) in &convs{if k>0&&k<=max_denom{best=(h,k);}}best
}

/// Continued fraction for √n (periodic part).
pub fn sqrt_cf(n:u64)->(u64,Vec<u64>){
    let a0=crate::powers::sqrt(n as f64)as u64;if a0*a0==n{return(a0,vec![]);}
    let mut period=Vec::new();let(mut m,mut d,mut a)=(0u64,1u64,a0);
    loop{m=d*a-m;d=(n-m*m)/d;a=(a0+m)/d;period.push(a);if a==2*a0{break;}}
    (a0,period)
}

/// Solve Pell's equation x² - n*y² = 1 using continued fractions.
pub fn pell_fundamental(n:u64)->Option<(u64,u64)>{
    let(a0,period)=sqrt_cf(n);if period.is_empty(){return None;}
    let mut cf=vec![a0];for &p in &period{cf.push(p);}
    // Try convergents until we find a solution
    let len=period.len();
    let idx=if len%2==0{len-1}else{2*len-1};
    let mut full_cf=vec![a0];
    for i in 0..=idx{full_cf.push(period[i%len]);}
    let convs=convergents(&full_cf.iter().map(|&x|x as i64).collect::<Vec<_>>());
    for &(h,k) in &convs{if h>0&&k>0{let(x,y)=(h as u64,k as u64);
        if x*x==n*y*y+1{return Some((x,y));}}}None
}

/// Golden ratio via continued fraction [1; 1, 1, 1, ...].
pub fn golden_ratio_cf(terms:usize)->f64{evaluate_cf(&vec![1i64;terms])}

/// e via continued fraction [2; 1, 2, 1, 1, 4, 1, 1, 6, ...].
pub fn e_cf(terms:usize)->f64{
    let mut cf=vec![2i64];let mut k=1;
    while cf.len()<terms{cf.push(1);cf.push(2*k);cf.push(1);k+=1;}
    cf.truncate(terms);evaluate_cf(&cf)
}

/// π via continued fraction (Stern's).
pub fn pi_cf_approx()->(i64,i64){best_rational(crate::constants::PI,10000)}

/// Generalized continued fraction: a0 + b1/(a1 + b2/(a2 + ...)).
pub fn generalized_cf(a:&[f64],b:&[f64])->f64{
    let n=a.len();let mut result=a[n-1];
    for i in (0..n-1).rev(){result=a[i]+b[i+1]/result;}result
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_pi(){let cf=continued_fraction(crate::constants::PI,10);assert_eq!(cf[0],3);assert_eq!(cf[1],7);}
    #[test]fn test_convergent(){let cf=continued_fraction(crate::constants::PI,4);let c=convergents(&cf);
        let(h,k)=c[1];assert_eq!((h,k),(22,7));} // 22/7
    #[test]fn test_golden(){assert!(approx_eq(golden_ratio_cf(30),1.618033988,1e-6));}
    #[test]fn test_sqrt2_cf(){let(a0,period)=sqrt_cf(2);assert_eq!(a0,1);assert_eq!(period,vec![2]);}
    #[test]fn test_pell(){let(x,y)=pell_fundamental(2).unwrap();assert_eq!((x,y),(3,2));}
}
