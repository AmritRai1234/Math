//! Bézier curves, parametric curves, arc length.

use crate::powers;

/// Evaluate cubic Bézier curve at parameter t ∈ [0,1].
pub fn bezier_cubic(p0:(f64,f64),p1:(f64,f64),p2:(f64,f64),p3:(f64,f64),t:f64)->(f64,f64){
    let u=1.0-t;let(u2,u3)=(u*u,u*u*u);let(t2,t3)=(t*t,t*t*t);
    (u3*p0.0+3.0*u2*t*p1.0+3.0*u*t2*p2.0+t3*p3.0,u3*p0.1+3.0*u2*t*p1.1+3.0*u*t2*p2.1+t3*p3.1)
}

/// Evaluate quadratic Bézier.
pub fn bezier_quadratic(p0:(f64,f64),p1:(f64,f64),p2:(f64,f64),t:f64)->(f64,f64){
    let u=1.0-t;(u*u*p0.0+2.0*u*t*p1.0+t*t*p2.0,u*u*p0.1+2.0*u*t*p1.1+t*t*p2.1)
}

/// General Bézier curve (de Casteljau algorithm).
pub fn bezier(points:&[(f64,f64)],t:f64)->(f64,f64){
    let mut pts:Vec<(f64,f64)>=points.to_vec();
    while pts.len()>1{pts=(0..pts.len()-1).map(|i|{let u=1.0-t;
        (u*pts[i].0+t*pts[i+1].0,u*pts[i].1+t*pts[i+1].1)}).collect();}pts[0]
}

/// Bézier curve derivative at t.
pub fn bezier_derivative(points:&[(f64,f64)],t:f64)->(f64,f64){
    let n=points.len()-1;let dpts:Vec<(f64,f64)>=(0..n).map(|i|
        (n as f64*(points[i+1].0-points[i].0),n as f64*(points[i+1].1-points[i].1))).collect();
    bezier(&dpts,t)
}

/// Arc length of parametric curve using numerical integration.
pub fn arc_length(fx:&dyn Fn(f64)->f64,fy:&dyn Fn(f64)->f64,a:f64,b:f64,n:usize)->f64{
    let dt=(b-a)/n as f64;let mut length=0.0;
    for i in 0..n{let t0=a+i as f64*dt;let t1=t0+dt;
        let(dx,dy)=(fx(t1)-fx(t0),fy(t1)-fy(t0));length+=powers::sqrt(dx*dx+dy*dy);}length
}

/// Arc length of Bézier curve (numerical).
pub fn bezier_arc_length(points:&[(f64,f64)],n:usize)->f64{
    let mut length=0.0;let mut prev=points[0];
    for i in 1..=n{let t=i as f64/n as f64;let p=bezier(points,t);
        let(dx,dy)=(p.0-prev.0,p.1-prev.1);length+=powers::sqrt(dx*dx+dy*dy);prev=p;}length
}

/// Curvature of parametric curve at point (x'(t), y'(t), x''(t), y''(t)).
pub fn parametric_curvature(dx:f64,dy:f64,ddx:f64,ddy:f64)->f64{
    let num=crate::arithmetic::abs(dx*ddy-dy*ddx);let den=powers::pow(dx*dx+dy*dy,1.5);
    if den<1e-30{0.0}else{num/den}
}

/// Subdivide Bézier curve at t (de Casteljau).
pub fn bezier_subdivide(points:&[(f64,f64)],t:f64)->(Vec<(f64,f64)>,Vec<(f64,f64)>){
    let n=points.len();let mut left=Vec::new();let mut right=Vec::new();
    let mut levels=vec![points.to_vec()];
    for _ in 1..n{let prev=levels.last().unwrap();
        let next:Vec<(f64,f64)>=(0..prev.len()-1).map(|i|{let u=1.0-t;
            (u*prev[i].0+t*prev[i+1].0,u*prev[i].1+t*prev[i+1].1)}).collect();
        levels.push(next);}
    for level in &levels{left.push(level[0]);}
    for level in levels.iter().rev(){right.push(*level.last().unwrap());}
    (left,right)
}

/// Catmull-Rom spline through points.
pub fn catmull_rom(p0:(f64,f64),p1:(f64,f64),p2:(f64,f64),p3:(f64,f64),t:f64)->(f64,f64){
    let t2=t*t;let t3=t2*t;
    let x=0.5*((2.0*p1.0)+(-p0.0+p2.0)*t+(2.0*p0.0-5.0*p1.0+4.0*p2.0-p3.0)*t2+(-p0.0+3.0*p1.0-3.0*p2.0+p3.0)*t3);
    let y=0.5*((2.0*p1.1)+(-p0.1+p2.1)*t+(2.0*p0.1-5.0*p1.1+4.0*p2.1-p3.1)*t2+(-p0.1+3.0*p1.1-3.0*p2.1+p3.1)*t3);
    (x,y)
}

/// Lissajous curve point.
pub fn lissajous(a:f64,b:f64,delta:f64,t:f64)->(f64,f64){
    (crate::trig::sin(a*t+delta),crate::trig::sin(b*t))
}

/// Parametric circle.
pub fn circle(r:f64,t:f64)->(f64,f64){(r*crate::trig::cos(t),r*crate::trig::sin(t))}

/// Parametric ellipse.
pub fn ellipse(a:f64,b:f64,t:f64)->(f64,f64){(a*crate::trig::cos(t),b*crate::trig::sin(t))}

/// Parametric spiral (Archimedean).
pub fn archimedean_spiral(a:f64,b:f64,t:f64)->(f64,f64){
    let r=a+b*t;(r*crate::trig::cos(t),r*crate::trig::sin(t))
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_bezier_endpoints(){let p=bezier_cubic((0.0,0.0),(1.0,1.0),(2.0,1.0),(3.0,0.0),0.0);
        assert!(approx_eq(p.0,0.0,1e-10));assert!(approx_eq(p.1,0.0,1e-10));}
    #[test]fn test_bezier_end(){let p=bezier_cubic((0.0,0.0),(1.0,1.0),(2.0,1.0),(3.0,0.0),1.0);
        assert!(approx_eq(p.0,3.0,1e-10));}
    #[test]fn test_decasteljau(){let p=bezier(&[(0.0,0.0),(1.0,1.0),(2.0,0.0)],0.5);
        assert!(approx_eq(p.0,1.0,1e-10));}
    #[test]fn test_arc_length_circle(){let l=arc_length(&|t|crate::trig::cos(t),&|t|crate::trig::sin(t),0.0,2.0*crate::constants::PI,10000);
        assert!(approx_eq(l,2.0*crate::constants::PI,0.01));}
}
