//! Complex number arithmetic and functions.

use crate::exponential;
use crate::powers;
use crate::trig;
use std::ops::{Add, Sub, Mul, Div, Neg};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Complex { pub re: f64, pub im: f64 }

impl Complex {
    pub fn new(re:f64,im:f64)->Self{Self{re,im}}
    pub fn from_polar(r:f64,theta:f64)->Self{Self{re:r*trig::cos(theta),im:r*trig::sin(theta)}}
    pub fn magnitude(self)->f64{powers::sqrt(self.re*self.re+self.im*self.im)}
    pub fn magnitude_sq(self)->f64{self.re*self.re+self.im*self.im}
    pub fn phase(self)->f64{trig::atan2(self.im,self.re)}
    pub fn to_polar(self)->(f64,f64){(self.magnitude(),self.phase())}
    pub fn conjugate(self)->Self{Self{re:self.re,im:-self.im}}
    pub fn reciprocal(self)->Self{let m2=self.magnitude_sq();Self{re:self.re/m2,im:-self.im/m2}}

    pub fn exp(self)->Self{let e=exponential::exp(self.re);Self{re:e*trig::cos(self.im),im:e*trig::sin(self.im)}}
    pub fn ln(self)->Self{Self{re:exponential::ln(self.magnitude()),im:self.phase()}}
    pub fn pow(self,w:Self)->Self{if self.re==0.0&&self.im==0.0{return Self::new(0.0,0.0);}(self.ln()*w).exp()}
    pub fn sqrt(self)->Self{let r=powers::sqrt(self.magnitude());let theta=self.phase()/2.0;Self::from_polar(r,theta)}

    pub fn sin(self)->Self{Self{re:trig::sin(self.re)*crate::hyperbolic::cosh(self.im),im:trig::cos(self.re)*crate::hyperbolic::sinh(self.im)}}
    pub fn cos(self)->Self{Self{re:trig::cos(self.re)*crate::hyperbolic::cosh(self.im),im:-trig::sin(self.re)*crate::hyperbolic::sinh(self.im)}}

    pub fn nth_roots(self,n:u32)->Vec<Self>{
        let(r,theta)=self.to_polar();let r_n=powers::nth_root(r,n);
        (0..n).map(|k|{let angle=(theta+2.0*crate::constants::PI*k as f64)/n as f64;Self::from_polar(r_n,angle)}).collect()
    }
}

impl Add for Complex{type Output=Self;fn add(self,o:Self)->Self{Self{re:self.re+o.re,im:self.im+o.im}}}
impl Sub for Complex{type Output=Self;fn sub(self,o:Self)->Self{Self{re:self.re-o.re,im:self.im-o.im}}}
impl Mul for Complex{type Output=Self;fn mul(self,o:Self)->Self{Self{re:self.re*o.re-self.im*o.im,im:self.re*o.im+self.im*o.re}}}
impl Div for Complex{type Output=Self;fn div(self,o:Self)->Self{let d=o.magnitude_sq();Self{re:(self.re*o.re+self.im*o.im)/d,im:(self.im*o.re-self.re*o.im)/d}}}
impl Neg for Complex{type Output=Self;fn neg(self)->Self{Self{re:-self.re,im:-self.im}}}

impl std::fmt::Display for Complex{
    fn fmt(&self,f:&mut std::fmt::Formatter<'_>)->std::fmt::Result{
        if self.im>=0.0{write!(f,"{}+{}i",self.re,self.im)}else{write!(f,"{}{}i",self.re,self.im)}
    }
}

#[cfg(test)]
mod tests{
    use super::*;use crate::arithmetic::approx_eq;use crate::constants;
    #[test]fn test_arithmetic(){let a=Complex::new(1.0,2.0);let b=Complex::new(3.0,4.0);let c=a*b;assert!(approx_eq(c.re,-5.0,1e-12));assert!(approx_eq(c.im,10.0,1e-12));}
    #[test]fn test_polar(){let z=Complex::new(0.0,1.0);let(r,t)=z.to_polar();assert!(approx_eq(r,1.0,1e-12));assert!(approx_eq(t,constants::FRAC_PI_2,1e-12));}
    #[test]fn test_euler(){let z=Complex::new(0.0,constants::PI).exp();assert!(approx_eq(z.re,-1.0,1e-12));assert!(approx_eq(z.im,0.0,1e-12));}
    #[test]fn test_nth_roots(){let z=Complex::new(1.0,0.0);let roots=z.nth_roots(4);assert_eq!(roots.len(),4);}
}
