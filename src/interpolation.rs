//! Interpolation functions.

use crate::arithmetic::clamp;

/// Linear interpolation.
pub fn lerp(a:f64,b:f64,t:f64)->f64{a+t*(b-a)}

/// Bilinear interpolation.
pub fn bilerp(c00:f64,c10:f64,c01:f64,c11:f64,tx:f64,ty:f64)->f64{
    let a=lerp(c00,c10,tx);let b=lerp(c01,c11,tx);lerp(a,b,ty)
}

/// Cubic interpolation (Catmull-Rom style).
pub fn cubic_interp(p0:f64,p1:f64,p2:f64,p3:f64,t:f64)->f64{
    let t2=t*t;let t3=t2*t;
    0.5*((-p0+3.0*p1-3.0*p2+p3)*t3+(2.0*p0-5.0*p1+4.0*p2-p3)*t2+(-p0+p2)*t+2.0*p1)
}

/// Catmull-Rom spline (same as cubic_interp).
pub fn catmull_rom(p0:f64,p1:f64,p2:f64,p3:f64,t:f64)->f64{cubic_interp(p0,p1,p2,p3,t)}

/// Quadratic Bézier: B(t) = (1-t)²P0 + 2(1-t)tP1 + t²P2.
pub fn bezier2(p0:f64,p1:f64,p2:f64,t:f64)->f64{
    let u=1.0-t;u*u*p0+2.0*u*t*p1+t*t*p2
}

/// Cubic Bézier.
pub fn bezier3(p0:f64,p1:f64,p2:f64,p3:f64,t:f64)->f64{
    let u=1.0-t;u*u*u*p0+3.0*u*u*t*p1+3.0*u*t*t*p2+t*t*t*p3
}

/// Hermite interpolation.
pub fn hermite_interp(p0:f64,m0:f64,p1:f64,m1:f64,t:f64)->f64{
    let t2=t*t;let t3=t2*t;
    (2.0*t3-3.0*t2+1.0)*p0+(t3-2.0*t2+t)*m0+(-2.0*t3+3.0*t2)*p1+(t3-t2)*m1
}

/// Lagrange interpolation at x given points.
pub fn lagrange(points:&[(f64,f64)],x:f64)->f64{
    let n=points.len();let mut result=0.0;
    for i in 0..n{let mut basis=1.0;for j in 0..n{if i!=j{basis*=(x-points[j].0)/(points[i].0-points[j].0);}}result+=basis*points[i].1;}
    result
}

/// Newton's interpolation.
pub fn newton_interp(points:&[(f64,f64)],x:f64)->f64{
    let n=points.len();let mut dd=vec![0.0;n];for i in 0..n{dd[i]=points[i].1;}
    for j in 1..n{for i in (j..n).rev(){dd[i]=(dd[i]-dd[i-1])/(points[i].0-points[i-j].0);}}
    let mut result=dd[n-1];for i in (0..n-1).rev(){result=result*(x-points[i].0)+dd[i];}
    result
}

/// Smoothstep.
pub fn smoothstep(e0:f64,e1:f64,x:f64)->f64{let t=clamp((x-e0)/(e1-e0),0.0,1.0);t*t*(3.0-2.0*t)}

/// Smootherstep (Ken Perlin).
pub fn smootherstep(e0:f64,e1:f64,x:f64)->f64{let t=clamp((x-e0)/(e1-e0),0.0,1.0);t*t*t*(t*(t*6.0-15.0)+10.0)}

/// Inverse lerp.
pub fn inverse_lerp(a:f64,b:f64,v:f64)->f64{if(b-a).abs()<f64::EPSILON{0.0}else{(v-a)/(b-a)}}

#[cfg(test)]
mod tests{
    use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_lerp(){assert!(approx_eq(lerp(0.0,10.0,0.5),5.0,1e-14));}
    #[test]fn test_bezier(){assert!(approx_eq(bezier2(0.0,5.0,10.0,0.5),5.0,1e-14));}
    #[test]fn test_lagrange(){let pts=[(0.0,0.0),(1.0,1.0),(2.0,4.0)];assert!(approx_eq(lagrange(&pts,1.5),2.25,1e-12));}
}
