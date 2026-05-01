//! Exact rational number arithmetic.

use crate::number_theory::gcd;
use std::ops::{Add,Sub,Mul,Div,Neg};
use std::cmp::Ordering;

#[derive(Clone,Copy,Debug)]
pub struct Rational{pub num:i64,pub den:i64}

impl Rational{
    pub fn new(num:i64,den:i64)->Self{
        if den==0{panic!("Rational: zero denominator");}
        let sign=if den<0{-1}else{1};
        let g=gcd(num.unsigned_abs(),den.unsigned_abs())as i64;
        Self{num:sign*num/g,den:sign*den/g}
    }
    pub fn from_int(n:i64)->Self{Self{num:n,den:1}}
    pub fn zero()->Self{Self{num:0,den:1}}
    pub fn one()->Self{Self{num:1,den:1}}
    pub fn is_zero(self)->bool{self.num==0}
    pub fn is_integer(self)->bool{self.den==1}
    pub fn abs(self)->Self{Self{num:self.num.abs(),den:self.den}}
    pub fn recip(self)->Self{Self::new(self.den,self.num)}
    pub fn pow(self,n:u32)->Self{Self::new(self.num.pow(n),self.den.pow(n))}
    pub fn to_f64(self)->f64{self.num as f64/self.den as f64}
    pub fn mediant(self,other:Self)->Self{Self::new(self.num+other.num,self.den+other.den)}
    pub fn floor(self)->i64{if self.num>=0{self.num/self.den}else{(self.num-self.den+1)/self.den}}
    pub fn ceil(self)->i64{if self.num>=0{(self.num+self.den-1)/self.den}else{self.num/self.den}}

    /// Convert f64 to rational approximation using continued fractions.
    pub fn from_f64(x:f64,max_den:i64)->Self{
        let sign=if x<0.0{-1}else{1};let x=x.abs();
        let mut p0=0i64;let mut q0=1i64;let mut p1=1i64;let mut q1=0i64;
        let mut rem=x;
        for _ in 0..64{let a=rem as i64;
            let p2=a*p1+p0;let q2=a*q1+q0;
            if q2>max_den{break;}
            p0=p1;q0=q1;p1=p2;q1=q2;
            let frac=rem-a as f64;if frac.abs()<1e-12{break;}
            rem=1.0/frac;
        }
        Self::new(sign*p1,q1)
    }
}

impl PartialEq for Rational{fn eq(&self,o:&Self)->bool{self.num==o.num&&self.den==o.den}}
impl Eq for Rational{}
impl PartialOrd for Rational{fn partial_cmp(&self,o:&Self)->Option<Ordering>{Some(self.cmp(o))}}
impl Ord for Rational{fn cmp(&self,o:&Self)->Ordering{(self.num*o.den).cmp(&(o.num*self.den))}}

impl Add for Rational{type Output=Self;fn add(self,o:Self)->Self{Self::new(self.num*o.den+o.num*self.den,self.den*o.den)}}
impl Sub for Rational{type Output=Self;fn sub(self,o:Self)->Self{Self::new(self.num*o.den-o.num*self.den,self.den*o.den)}}
impl Mul for Rational{type Output=Self;fn mul(self,o:Self)->Self{Self::new(self.num*o.num,self.den*o.den)}}
impl Div for Rational{type Output=Self;fn div(self,o:Self)->Self{Self::new(self.num*o.den,self.den*o.num)}}
impl Neg for Rational{type Output=Self;fn neg(self)->Self{Self{num:-self.num,den:self.den}}}

impl std::fmt::Display for Rational{
    fn fmt(&self,f:&mut std::fmt::Formatter<'_>)->std::fmt::Result{
        if self.den==1{write!(f,"{}",self.num)}else{write!(f,"{}/{}",self.num,self.den)}
    }
}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_basic(){let a=Rational::new(1,2);let b=Rational::new(1,3);let c=a+b;assert_eq!(c,Rational::new(5,6));}
    #[test]fn test_simplify(){assert_eq!(Rational::new(4,6),Rational::new(2,3));}
    #[test]fn test_mul(){let r=Rational::new(2,3)*Rational::new(3,4);assert_eq!(r,Rational::new(1,2));}
    #[test]fn test_from_f64(){let r=Rational::from_f64(0.333333,1000);assert_eq!(r,Rational::new(1,3));}
    #[test]fn test_to_f64(){assert!((Rational::new(1,3).to_f64()-0.33333333).abs()<0.001);}
    #[test]fn test_ord(){assert!(Rational::new(1,3)<Rational::new(1,2));}
}
