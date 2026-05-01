//! Interval arithmetic with guaranteed bounds.

use std::ops::{Add,Sub,Mul,Div,Neg};

#[derive(Clone,Copy,Debug,PartialEq)]
pub struct Interval{pub lo:f64,pub hi:f64}

impl Interval{
    pub fn new(lo:f64,hi:f64)->Self{Self{lo:lo.min(hi),hi:lo.max(hi)}}
    pub fn point(x:f64)->Self{Self{lo:x,hi:x}}
    pub fn width(self)->f64{self.hi-self.lo}
    pub fn midpoint(self)->f64{(self.lo+self.hi)/2.0}
    pub fn radius(self)->f64{self.width()/2.0}
    pub fn contains(self,x:f64)->bool{x>=self.lo&&x<=self.hi}
    pub fn overlaps(self,o:Self)->bool{self.lo<=o.hi&&o.lo<=self.hi}
    pub fn hull(self,o:Self)->Self{Self{lo:self.lo.min(o.lo),hi:self.hi.max(o.hi)}}
    pub fn intersect(self,o:Self)->Option<Self>{if self.overlaps(o){Some(Self{lo:self.lo.max(o.lo),hi:self.hi.min(o.hi)})}else{None}}
    pub fn abs(self)->Self{if self.lo>=0.0{self}else if self.hi<=0.0{Self{lo:-self.hi,hi:-self.lo}}else{Self{lo:0.0,hi:self.lo.abs().max(self.hi.abs())}}}
    pub fn pow2(self)->Self{if self.lo>=0.0{Self{lo:self.lo*self.lo,hi:self.hi*self.hi}}else if self.hi<=0.0{Self{lo:self.hi*self.hi,hi:self.lo*self.lo}}else{Self{lo:0.0,hi:self.lo.abs().max(self.hi.abs()).powi(2)}}}
    pub fn sqrt(self)->Self{if self.lo<0.0{Self{lo:0.0,hi:crate::powers::sqrt(self.hi.max(0.0))}}else{Self{lo:crate::powers::sqrt(self.lo),hi:crate::powers::sqrt(self.hi)}}}
    pub fn exp(self)->Self{Self{lo:crate::exponential::exp(self.lo),hi:crate::exponential::exp(self.hi)}}
    pub fn ln(self)->Self{Self{lo:crate::exponential::ln(self.lo.max(f64::MIN_POSITIVE)),hi:crate::exponential::ln(self.hi.max(f64::MIN_POSITIVE))}}
    pub fn sin(self)->Self{// Conservative bound
        if self.width()>=crate::constants::TAU{return Self{lo:-1.0,hi:1.0};}
        let a=crate::trig::sin(self.lo);let b=crate::trig::sin(self.hi);
        let mut lo=a.min(b);let mut hi=a.max(b);
        // Check if interval contains critical points
        let check=|angle:f64|{let n_lo=((self.lo-angle)/crate::constants::TAU).floor();let n_hi=((self.hi-angle)/crate::constants::TAU).floor();n_lo!=n_hi};
        if check(crate::constants::FRAC_PI_2){hi=1.0;}
        if check(-crate::constants::FRAC_PI_2){lo=-1.0;}
        Self{lo,hi}
    }
    pub fn cos(self)->Self{Self{lo:self.lo+crate::constants::FRAC_PI_2,hi:self.hi+crate::constants::FRAC_PI_2}.sin()}

    /// Check if 0 is in the interval.
    pub fn contains_zero(self)->bool{self.lo<=0.0&&self.hi>=0.0}
}

impl Add for Interval{type Output=Self;fn add(self,o:Self)->Self{Self{lo:self.lo+o.lo,hi:self.hi+o.hi}}}
impl Sub for Interval{type Output=Self;fn sub(self,o:Self)->Self{Self{lo:self.lo-o.hi,hi:self.hi-o.lo}}}
impl Mul for Interval{type Output=Self;fn mul(self,o:Self)->Self{
    let products=[self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi];
    Self{lo:products.iter().cloned().fold(f64::INFINITY,f64::min),hi:products.iter().cloned().fold(f64::NEG_INFINITY,f64::max)}}}
impl Div for Interval{type Output=Self;fn div(self,o:Self)->Self{
    if o.contains_zero(){return Self{lo:f64::NEG_INFINITY,hi:f64::INFINITY};}
    self*Self{lo:1.0/o.hi,hi:1.0/o.lo}}}
impl Neg for Interval{type Output=Self;fn neg(self)->Self{Self{lo:-self.hi,hi:-self.lo}}}

impl std::fmt::Display for Interval{fn fmt(&self,f:&mut std::fmt::Formatter)->std::fmt::Result{write!(f,"[{}, {}]",self.lo,self.hi)}}

/// Interval Newton's method for root finding with guaranteed inclusion.
pub fn interval_newton(f:&dyn Fn(Interval)->Interval,df:&dyn Fn(Interval)->Interval,mut x:Interval,iters:usize)->Interval{
    for _ in 0..iters{let fx=f(Interval::point(x.midpoint()));let dfx=df(x);
        if dfx.contains_zero(){break;}
        let new=Interval::point(x.midpoint())-fx/dfx;
        if let Some(intersection)=x.intersect(new){x=intersection;}else{break;}}
    x
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_add(){let a=Interval::new(1.0,2.0);let b=Interval::new(3.0,4.0);let c=a+b;assert!(approx_eq(c.lo,4.0,1e-10));assert!(approx_eq(c.hi,6.0,1e-10));}
    #[test]fn test_mul(){let a=Interval::new(-1.0,2.0);let b=Interval::new(3.0,4.0);let c=a*b;assert!(approx_eq(c.lo,-4.0,1e-10));assert!(approx_eq(c.hi,8.0,1e-10));}
    #[test]fn test_contains(){assert!(Interval::new(1.0,3.0).contains(2.0));assert!(!Interval::new(1.0,3.0).contains(4.0));}
    #[test]fn test_sqrt(){let i=Interval::new(4.0,9.0).sqrt();assert!(approx_eq(i.lo,2.0,1e-10));assert!(approx_eq(i.hi,3.0,1e-10));}
    #[test]fn test_width(){assert!(approx_eq(Interval::new(1.0,5.0).width(),4.0,1e-10));}
}
