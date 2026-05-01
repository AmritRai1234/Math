//! Root finding and optimization.

use crate::arithmetic::abs;

pub fn bisection(f:&dyn Fn(f64)->f64,mut a:f64,mut b:f64,tol:f64)->f64{
    for _ in 0..1000{let mid=(a+b)/2.0;if abs(f(mid))<tol||abs(b-a)<tol{return mid;}if f(a)*f(mid)<0.0{b=mid;}else{a=mid;}}(a+b)/2.0
}

pub fn newton_raphson(f:&dyn Fn(f64)->f64,df:&dyn Fn(f64)->f64,mut x:f64,tol:f64)->f64{
    for _ in 0..1000{let fx=f(x);if abs(fx)<tol{return x;}let dfx=df(x);if abs(dfx)<1e-30{break;}x-=fx/dfx;}x
}

pub fn secant(f:&dyn Fn(f64)->f64,mut x0:f64,mut x1:f64,tol:f64)->f64{
    for _ in 0..1000{let f0=f(x0);let f1=f(x1);if abs(f1)<tol{return x1;}let d=f1-f0;if abs(d)<1e-30{break;}let x2=x1-f1*(x1-x0)/d;x0=x1;x1=x2;}x1
}

pub fn brent(f:&dyn Fn(f64)->f64,mut a:f64,mut b:f64,tol:f64)->f64{
    let mut fa=f(a);let mut fb=f(b);if fa*fb>0.0{return f64::NAN;}
    if abs(fa)<abs(fb){std::mem::swap(&mut a,&mut b);std::mem::swap(&mut fa,&mut fb);}
    let mut c=a;let mut _fc=fa;let mut d=b-a;
    for _ in 0..1000{if abs(fb)<tol{return b;}d=(a-b)/2.0;b+=d;fb=f(b);
        if fa*fb>0.0{a=c;fa=_fc;}c=b;_fc=fb;if abs(fa)<abs(fb){std::mem::swap(&mut a,&mut b);std::mem::swap(&mut fa,&mut fb);}}b
}

pub fn fixed_point(g:&dyn Fn(f64)->f64,mut x:f64,tol:f64)->f64{
    for _ in 0..10000{let next=g(x);if abs(next-x)<tol{return next;}x=next;}x
}

pub fn gradient_descent_1d(f:&dyn Fn(f64)->f64,mut x:f64,lr:f64,steps:usize)->f64{
    let h=1e-8;for _ in 0..steps{let grad=(f(x+h)-f(x-h))/(2.0*h);x-=lr*grad;}x
}

pub fn golden_section(f:&dyn Fn(f64)->f64,mut a:f64,mut b:f64,tol:f64)->f64{
    let gr=(crate::powers::sqrt(5.0)-1.0)/2.0;
    let mut c=b-gr*(b-a);let mut d=a+gr*(b-a);
    for _ in 0..1000{if abs(b-a)<tol{break;}if f(c)<f(d){b=d;}else{a=c;}c=b-gr*(b-a);d=a+gr*(b-a);}(a+b)/2.0
}

pub fn nelder_mead_2d(f:&dyn Fn(f64,f64)->f64,start:(f64,f64),step:f64,tol:f64)->(f64,f64){
    let mut s=[(start.0,start.1),(start.0+step,start.1),(start.0,start.1+step)];
    let mut fv=[f(s[0].0,s[0].1),f(s[1].0,s[1].1),f(s[2].0,s[2].1)];
    for _ in 0..10000{let mut idx=[0,1,2];idx.sort_by(|&a,&b|fv[a].partial_cmp(&fv[b]).unwrap());
        if abs(fv[idx[2]]-fv[idx[0]])<tol{return s[idx[0]];}
        let cx=(s[idx[0]].0+s[idx[1]].0)/2.0;let cy=(s[idx[0]].1+s[idx[1]].1)/2.0;
        let rx=2.0*cx-s[idx[2]].0;let ry=2.0*cy-s[idx[2]].1;let fr=f(rx,ry);
        if fr<fv[idx[1]]{s[idx[2]]=(rx,ry);fv[idx[2]]=fr;}
        else{for i in 1..3{s[idx[i]]=((s[idx[0]].0+s[idx[i]].0)/2.0,(s[idx[0]].1+s[idx[i]].1)/2.0);fv[idx[i]]=f(s[idx[i]].0,s[idx[i]].1);}}
    }s[0]
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_bisection(){let r=bisection(&|x|x*x-2.0,1.0,2.0,1e-12);assert!(approx_eq(r,crate::powers::sqrt(2.0),1e-10));}
    #[test]fn test_newton(){let r=newton_raphson(&|x|x*x-2.0,&|x|2.0*x,1.5,1e-12);assert!(approx_eq(r,crate::powers::sqrt(2.0),1e-10));}
    #[test]fn test_golden(){let r=golden_section(&|x|(x-3.0)*(x-3.0),0.0,10.0,1e-8);assert!(approx_eq(r,3.0,1e-5));}
}
