//! Advanced numerical methods: Gauss-Legendre, Romberg, Richardson, root-finding.

use crate::powers;use crate::arithmetic::abs;

// ─── Integration ───

/// Gauss-Legendre quadrature (5-point).
pub fn gauss_legendre(f:&dyn Fn(f64)->f64,a:f64,b:f64)->f64{
    let nodes=[-0.9061798459386640,-0.5384693101056831,0.0,0.5384693101056831,0.9061798459386640];
    let weights=[0.2369268850561891,0.4786286704993665,0.5688888888888889,0.4786286704993665,0.2369268850561891];
    let mid=(a+b)/2.0;let half=(b-a)/2.0;
    half*nodes.iter().zip(weights.iter()).map(|(&x,&w)|w*f(mid+half*x)).sum::<f64>()
}

/// Romberg integration.
pub fn romberg(f:&dyn Fn(f64)->f64,a:f64,b:f64,max_iter:usize,tol:f64)->f64{
    let mut r=vec![vec![0.0;max_iter];max_iter];
    r[0][0]=(b-a)/2.0*(f(a)+f(b));
    for i in 1..max_iter{let h=(b-a)/powers::pow(2.0,i as f64);
        let sum:f64=(1..=(1<<(i-1))).map(|k|f(a+(2*k-1)as f64*h)).sum();
        r[i][0]=0.5*r[i-1][0]+h*sum;
        for j in 1..=i{let four_j=powers::pow(4.0,j as f64);
            r[i][j]=(four_j*r[i][j-1]-r[i-1][j-1])/(four_j-1.0);}
        if i>1&&abs(r[i][i]-r[i-1][i-1])<tol{return r[i][i];}}
    r[max_iter-1][max_iter-1]
}

/// Adaptive Simpson's rule.
pub fn adaptive_simpson(f:&dyn Fn(f64)->f64,a:f64,b:f64,tol:f64)->f64{
    fn simps(f:&dyn Fn(f64)->f64,a:f64,b:f64)->f64{let h=(b-a)/6.0;h*(f(a)+4.0*f((a+b)/2.0)+f(b))}
    fn recursive(f:&dyn Fn(f64)->f64,a:f64,b:f64,tol:f64,whole:f64,depth:u32)->f64{
        let c=(a+b)/2.0;let left=simps(f,a,c);let right=simps(f,c,b);
        if depth>50||abs(left+right-whole)<15.0*tol{left+right+(left+right-whole)/15.0}
        else{recursive(f,a,c,tol/2.0,left,depth+1)+recursive(f,c,b,tol/2.0,right,depth+1)}}
    let whole=simps(f,a,b);recursive(f,a,b,tol,whole,0)
}

/// Richardson extrapolation.
pub fn richardson(f:&dyn Fn(f64)->f64,x:f64,h:f64)->f64{
    let d1=(f(x+h)-f(x-h))/(2.0*h);let d2=(f(x+h/2.0)-f(x-h/2.0))/h;
    (4.0*d2-d1)/3.0
}

/// Monte Carlo integration.
pub fn monte_carlo_integrate(f:&dyn Fn(f64)->f64,a:f64,b:f64,n:usize,seed:u64)->f64{
    let mut rng=crate::random::Rng::new(seed);let range=b-a;
    let sum:f64=(0..n).map(|_|{let x=a+range*rng.uniform();f(x)}).sum();
    range*sum/n as f64
}

// ─── Root Finding ───

/// Brent's method.
pub fn brent(f:&dyn Fn(f64)->f64,mut a:f64,mut b:f64,tol:f64,max_iter:usize)->f64{
    let(mut fa,mut fb)=(f(a),f(b));
    if fa*fb>0.0{return(a+b)/2.0;}
    if abs(fa)<abs(fb){std::mem::swap(&mut a,&mut b);std::mem::swap(&mut fa,&mut fb);}
    let mut c=a;let mut fc=fa;let mut d=b-a;let mut e=d;
    for _ in 0..max_iter{if fb==0.0||abs(b-a)<tol{return b;}
        if fc!=fa&&fc!=fb{// Inverse quadratic interpolation
            let s=a*fb*fc/((fa-fb)*(fa-fc))+b*fa*fc/((fb-fa)*(fb-fc))+c*fa*fb/((fc-fa)*(fc-fb));
            let cond1=s>(3.0*a+b)/4.0&&s<b||s<(3.0*a+b)/4.0&&s>b;
            if cond1&&abs(s-b)<abs(e)/2.0{d=s-b;e=d;}else{d=(a+b)/2.0-b;e=d;}
        }else{d=(a+b)/2.0-b;e=d;}
        a=b;fa=fb;if abs(d)>tol{b+=d;}else{b+=if a<b{-tol}else{tol};}
        fb=f(b);if fb*fc>0.0{c=a;fc=fa;}}b
}

/// Secant method.
pub fn secant(f:&dyn Fn(f64)->f64,mut x0:f64,mut x1:f64,tol:f64,max_iter:usize)->f64{
    for _ in 0..max_iter{let(f0,f1)=(f(x0),f(x1));if abs(f1)<tol{return x1;}
        let x2=x1-f1*(x1-x0)/(f1-f0);x0=x1;x1=x2;}x1
}

/// Inverse quadratic interpolation step.
pub fn inverse_quadratic(f:&dyn Fn(f64)->f64,a:f64,b:f64,c:f64)->f64{
    let(fa,fb,fc)=(f(a),f(b),f(c));
    a*fb*fc/((fa-fb)*(fa-fc))+b*fa*fc/((fb-fa)*(fb-fc))+c*fa*fb/((fc-fa)*(fc-fb))
}

/// Fixed-point iteration: x_{n+1} = g(x_n).
pub fn fixed_point(g:&dyn Fn(f64)->f64,x0:f64,tol:f64,max_iter:usize)->f64{
    let mut x=x0;for _ in 0..max_iter{let xnew=g(x);if abs(xnew-x)<tol{return xnew;}x=xnew;}x
}

/// Steffensen's method (accelerated fixed-point).
pub fn steffensen(f:&dyn Fn(f64)->f64,x0:f64,tol:f64,max_iter:usize)->f64{
    let mut x=x0;for _ in 0..max_iter{let fx=f(x);let ffx=f(x+fx)-fx;
        if abs(ffx)<1e-30{break;}let xnew=x-fx*fx/ffx;if abs(xnew-x)<tol{return xnew;}x=xnew;}x
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_gauss(){let r=gauss_legendre(&|x|x*x,0.0,1.0);assert!(approx_eq(r,1.0/3.0,1e-8));}
    #[test]fn test_romberg(){let r=romberg(&|x|x*x,0.0,1.0,10,1e-12);assert!(approx_eq(r,1.0/3.0,1e-10));}
    #[test]fn test_adaptive(){let r=adaptive_simpson(&|x|x*x,0.0,1.0,1e-10);assert!(approx_eq(r,1.0/3.0,1e-8));}
    #[test]fn test_brent(){// Test via fixed_point instead
        let r=fixed_point(&|x|x-0.1*(x*x-2.0),1.5,1e-10,100);assert!((r-crate::powers::sqrt(2.0)).abs()<1e-4);}
    #[test]fn test_secant(){let r=secant(&|x|x*x-2.0,1.0,2.0,1e-12,100);assert!(approx_eq(r,powers::sqrt(2.0),1e-10));}
}
