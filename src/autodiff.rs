//! Forward-mode automatic differentiation.

use std::ops::{Add,Sub,Mul,Div,Neg};

/// Dual number for forward-mode AD: f(a+bε) = f(a) + f'(a)bε.
#[derive(Clone,Copy,Debug,PartialEq)]
pub struct Dual{pub val:f64,pub der:f64}

impl Dual{
    pub fn new(val:f64,der:f64)->Self{Self{val,der}}
    pub fn constant(val:f64)->Self{Self{val,der:0.0}}
    pub fn variable(val:f64)->Self{Self{val,der:1.0}}

    pub fn sin(self)->Self{Self{val:crate::trig::sin(self.val),der:self.der*crate::trig::cos(self.val)}}
    pub fn cos(self)->Self{Self{val:crate::trig::cos(self.val),der:-self.der*crate::trig::sin(self.val)}}
    pub fn tan(self)->Self{let c=crate::trig::cos(self.val);Self{val:crate::trig::tan(self.val),der:self.der/(c*c)}}
    pub fn exp(self)->Self{let e=crate::exponential::exp(self.val);Self{val:e,der:self.der*e}}
    pub fn ln(self)->Self{Self{val:crate::exponential::ln(self.val),der:self.der/self.val}}
    pub fn sqrt(self)->Self{let s=crate::powers::sqrt(self.val);Self{val:s,der:self.der/(2.0*s)}}
    pub fn pow(self,n:f64)->Self{Self{val:crate::powers::pow(self.val,n),der:self.der*n*crate::powers::pow(self.val,n-1.0)}}
    pub fn abs(self)->Self{Self{val:crate::arithmetic::abs(self.val),der:if self.val>=0.0{self.der}else{-self.der}}}
    pub fn recip(self)->Self{Self{val:1.0/self.val,der:-self.der/(self.val*self.val)}}

    pub fn sinh(self)->Self{Self{val:crate::hyperbolic::sinh(self.val),der:self.der*crate::hyperbolic::cosh(self.val)}}
    pub fn cosh(self)->Self{Self{val:crate::hyperbolic::cosh(self.val),der:self.der*crate::hyperbolic::sinh(self.val)}}
    pub fn tanh(self)->Self{let t=crate::hyperbolic::tanh(self.val);Self{val:t,der:self.der*(1.0-t*t)}}

    pub fn asin(self)->Self{Self{val:crate::trig::asin(self.val),der:self.der/crate::powers::sqrt(1.0-self.val*self.val)}}
    pub fn acos(self)->Self{Self{val:crate::trig::acos(self.val),der:-self.der/crate::powers::sqrt(1.0-self.val*self.val)}}
    pub fn atan(self)->Self{Self{val:crate::trig::atan(self.val),der:self.der/(1.0+self.val*self.val)}}
}

impl Add for Dual{type Output=Self;fn add(self,o:Self)->Self{Self{val:self.val+o.val,der:self.der+o.der}}}
impl Sub for Dual{type Output=Self;fn sub(self,o:Self)->Self{Self{val:self.val-o.val,der:self.der-o.der}}}
impl Mul for Dual{type Output=Self;fn mul(self,o:Self)->Self{Self{val:self.val*o.val,der:self.val*o.der+self.der*o.val}}}
impl Div for Dual{type Output=Self;fn div(self,o:Self)->Self{Self{val:self.val/o.val,der:(self.der*o.val-self.val*o.der)/(o.val*o.val)}}}
impl Neg for Dual{type Output=Self;fn neg(self)->Self{Self{val:-self.val,der:-self.der}}}

impl Add<f64> for Dual{type Output=Self;fn add(self,o:f64)->Self{Self{val:self.val+o,der:self.der}}}
impl Mul<f64> for Dual{type Output=Self;fn mul(self,o:f64)->Self{Self{val:self.val*o,der:self.der*o}}}

/// Compute derivative of f at x using dual numbers.
pub fn derivative(f:&dyn Fn(Dual)->Dual,x:f64)->f64{f(Dual::variable(x)).der}

/// Compute gradient of f:Rⁿ→R using forward-mode AD.
pub fn gradient(f:&dyn Fn(&[Dual])->Dual,point:&[f64])->Vec<f64>{
    let n=point.len();(0..n).map(|i|{
        let args:Vec<Dual>=point.iter().enumerate().map(|(j,&v)|if j==i{Dual::variable(v)}else{Dual::constant(v)}).collect();
        f(&args).der
    }).collect()
}

/// Second derivative via nested dual (conceptual - uses finite diff on AD derivative).
pub fn second_derivative(f:&dyn Fn(Dual)->Dual,x:f64)->f64{
    let h=1e-6;(derivative(f,x+h)-derivative(f,x-h))/(2.0*h)
}

/// Jacobian of vector function f:Rⁿ→Rᵐ.
pub fn jacobian(fs:&[&dyn Fn(&[Dual])->Dual],point:&[f64])->Vec<Vec<f64>>{
    fs.iter().map(|fi|gradient(fi,point)).collect()
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_derivative_x2(){// d/dx(x²) = 2x at x=3 → 6
        assert!(approx_eq(derivative(&|x|x*x,3.0),6.0,1e-10));}
    #[test]fn test_derivative_sin(){// d/dx(sin(x)) = cos(x) at x=0 → 1
        assert!(approx_eq(derivative(&|x|x.sin(),0.0),1.0,1e-10));}
    #[test]fn test_chain_rule(){// d/dx(sin(x²)) = 2x*cos(x²) at x=1
        let expected=2.0*crate::trig::cos(1.0);
        assert!(approx_eq(derivative(&|x|(x*x).sin(),1.0),expected,1e-10));}
    #[test]fn test_gradient(){let f=|x:&[Dual]|x[0]*x[0]+x[1]*x[1];
        let g=gradient(&f,&[3.0,4.0]);assert!(approx_eq(g[0],6.0,1e-10));assert!(approx_eq(g[1],8.0,1e-10));}
    #[test]fn test_quotient(){// d/dx(x/(1+x²)) at x=1
        let d=derivative(&|x|x/(Dual::constant(1.0)+x*x),1.0);assert!(approx_eq(d,0.0,1e-10));}
}
