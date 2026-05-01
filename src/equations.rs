//! Algebraic equation solvers: quadratic, cubic, quartic, polynomial root-finding.

use crate::powers;use crate::trig;use crate::constants::PI;use crate::complex::Complex;

/// Quadratic formula: ax² + bx + c = 0. Returns (x1, x2) as complex.
pub fn quadratic(a:f64,b:f64,c:f64)->(Complex,Complex){
    let disc=b*b-4.0*a*c;
    if disc>=0.0{let sq=powers::sqrt(disc);
        (Complex::new((-b+sq)/(2.0*a),0.0),Complex::new((-b-sq)/(2.0*a),0.0))
    }else{let sq=powers::sqrt(-disc);
        (Complex::new(-b/(2.0*a),sq/(2.0*a)),Complex::new(-b/(2.0*a),-sq/(2.0*a)))}
}

/// Discriminant of quadratic.
pub fn quadratic_discriminant(a:f64,b:f64,c:f64)->f64{b*b-4.0*a*c}

/// Cubic equation: ax³ + bx² + cx + d = 0. Returns 3 real roots (Cardano's formula).
pub fn cubic(a:f64,b:f64,c:f64,d:f64)->Vec<f64>{
    let(b,c,d)=(b/a,c/a,d/a);
    let p=(3.0*c-b*b)/3.0;let q=(2.0*b*b*b-9.0*b*c+27.0*d)/27.0;
    let disc=q*q/4.0+p*p*p/27.0;
    if disc>0.0{// One real root
        let u=cbrt(-q/2.0+powers::sqrt(disc));let v=cbrt(-q/2.0-powers::sqrt(disc));
        vec![u+v-b/3.0]
    }else if disc.abs()<1e-14{// Repeated root
        let u=cbrt(-q/2.0);vec![2.0*u-b/3.0,-u-b/3.0]
    }else{// Three real roots (trigonometric method)
        let m=2.0*powers::sqrt(-p/3.0);let theta=trig::acos(3.0*q/(p*m))/3.0;
        vec![m*trig::cos(theta)-b/3.0,m*trig::cos(theta-2.0*PI/3.0)-b/3.0,m*trig::cos(theta-4.0*PI/3.0)-b/3.0]
    }
}

fn cbrt(x:f64)->f64{if x>=0.0{powers::pow(x,1.0/3.0)}else{-powers::pow(-x,1.0/3.0)}}

/// Quartic equation: ax⁴ + bx³ + cx² + dx + e = 0. Returns real roots via Ferrari's method.
pub fn quartic(a:f64,b:f64,c:f64,d:f64,e:f64)->Vec<f64>{
    let(b,c,d,e)=(b/a,c/a,d/a,e/a);
    // Depressed quartic: t⁴ + pt² + qt + r = 0
    let p=c-3.0*b*b/8.0;let q=b*b*b/8.0-b*c/2.0+d;
    let r=-3.0*b*b*b*b/256.0+b*b*c/16.0-b*d/4.0+e;
    // Resolvent cubic: y³ - p/2 y² - ry + (4pr - q²)/8 = 0
    let resolvent=cubic(1.0,-p/2.0,-r,(4.0*p*r-q*q)/8.0);
    let y=resolvent[0];
    let sq2y_p=2.0*y-p;if sq2y_p<0.0{return vec![];}
    let sq=powers::sqrt(sq2y_p);
    let mut roots=Vec::new();
    // Solve two quadratics
    for sign in [-1.0,1.0]{let disc=-2.0*y-p+sign*2.0*q/sq;
        if disc>=0.0{let s=powers::sqrt(disc);
            roots.push((-b/4.0+sign*sq/2.0+s/2.0).max(-b/4.0+sign*sq/2.0-s/2.0));
            roots.push((-b/4.0+sign*sq/2.0+s/2.0).min(-b/4.0+sign*sq/2.0-s/2.0));}}
    roots.sort_by(|a,b|a.partial_cmp(b).unwrap());roots.dedup();roots
}

/// Newton-Raphson for polynomial.
pub fn polynomial_newton(coeffs:&[f64],x0:f64,tol:f64,max_iter:usize)->f64{
    let eval=|coeffs:&[f64],x:f64|->f64{coeffs.iter().rev().fold(0.0,|acc,&c|acc*x+c)};
    let deriv_coeffs:Vec<f64>=coeffs.iter().enumerate().skip(1).map(|(i,&c)|c*i as f64).collect();
    let mut x=x0;for _ in 0..max_iter{let fx=eval(coeffs,x);let fpx=eval(&deriv_coeffs,x);
        if fpx.abs()<1e-30{break;}let dx=fx/fpx;x-=dx;if dx.abs()<tol{break;}}x
}

/// Durand-Kerner method: find all roots of polynomial simultaneously.
pub fn durand_kerner(coeffs:&[f64],max_iter:usize)->Vec<Complex>{
    let n=coeffs.len()-1;if n==0{return vec![];}
    let eval=|x:Complex|->Complex{coeffs.iter().enumerate().fold(Complex::new(0.0,0.0),|acc,(i,&c)|{
        let mut xk=Complex::new(1.0,0.0);for _ in 0..i{xk=Complex::new(xk.re*x.re-xk.im*x.im,xk.re*x.im+xk.im*x.re);}
        Complex::new(acc.re+c*xk.re,acc.im+c*xk.im)})};
    // Initial guesses on unit circle
    let mut roots:Vec<Complex>=(0..n).map(|k|{let angle=2.0*PI*k as f64/n as f64+0.4;
        Complex::new(0.7*trig::cos(angle),0.7*trig::sin(angle))}).collect();
    for _ in 0..max_iter{for i in 0..n{let fi=eval(roots[i]);
        let mut denom=Complex::new(coeffs[n],0.0);
        for j in 0..n{if j!=i{let d=Complex::new(roots[i].re-roots[j].re,roots[i].im-roots[j].im);
            denom=Complex::new(denom.re*d.re-denom.im*d.im,denom.re*d.im+denom.im*d.re);}}
        let mag=denom.re*denom.re+denom.im*denom.im;if mag>1e-30{
            roots[i]=Complex::new(roots[i].re-fi.re*denom.re/(mag)-fi.im*(-denom.im)/(mag),
                roots[i].im-(fi.im*denom.re-fi.re*(-denom.im))/(mag));}}}
    roots
}

/// Vieta's formulas: sum and product of roots from polynomial coefficients.
pub fn vieta_sum(coeffs:&[f64])->f64{let n=coeffs.len();-coeffs[n-2]/coeffs[n-1]}
pub fn vieta_product(coeffs:&[f64])->f64{let n=coeffs.len();let sign=if n%2==0{-1.0}else{1.0};sign*coeffs[0]/coeffs[n-1]}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_quadratic(){let(r1,r2)=quadratic(1.0,-3.0,2.0);// x²-3x+2=0 → x=1,2
        assert!(approx_eq(r1.re*r2.re,2.0,1e-10));}
    #[test]fn test_cubic(){let roots=cubic(1.0,-6.0,11.0,-6.0);// (x-1)(x-2)(x-3)=0
        assert!(roots.len()==3);}
    #[test]fn test_newton(){let r=polynomial_newton(&[-6.0,11.0,-6.0,1.0],0.5,1e-12,100);
        assert!(approx_eq(r,1.0,1e-6));}
}
