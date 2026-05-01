//! Spline interpolation: cubic splines, B-splines.

/// Natural cubic spline through points. Returns spline coefficients.
pub fn cubic_spline(x:&[f64],y:&[f64])->Vec<[f64;4]>{
    let n=x.len()-1;if n<1{return vec![];}
    let mut h=vec![0.0;n];for i in 0..n{h[i]=x[i+1]-x[i];}
    // Solve tridiagonal system for c coefficients
    let mut alpha=vec![0.0;n+1];
    for i in 1..n{alpha[i]=3.0*(y[i+1]-y[i])/h[i]-3.0*(y[i]-y[i-1])/h[i-1];}
    let mut l=vec![1.0;n+1];let mut mu=vec![0.0;n+1];let mut z=vec![0.0;n+1];
    for i in 1..n{l[i]=2.0*(x[i+1]-x[i-1])-h[i-1]*mu[i-1];
        mu[i]=h[i]/l[i];z[i]=(alpha[i]-h[i-1]*z[i-1])/l[i];}
    let mut c=vec![0.0;n+1];let mut b=vec![0.0;n];let mut d=vec![0.0;n];
    for j in(0..n).rev(){c[j]=z[j]-mu[j]*c[j+1];
        b[j]=(y[j+1]-y[j])/h[j]-h[j]*(c[j+1]+2.0*c[j])/3.0;
        d[j]=(c[j+1]-c[j])/(3.0*h[j]);}
    (0..n).map(|i|[y[i],b[i],c[i],d[i]]).collect()
}

/// Evaluate cubic spline at point t.
pub fn cubic_spline_eval(x:&[f64],coeffs:&[[f64;4]],t:f64)->f64{
    let n=coeffs.len();
    let mut i=0;for k in 0..n{if t>=x[k]{i=k;}}
    let dx=t-x[i];let[a,b,c,d]=coeffs[i];
    a+b*dx+c*dx*dx+d*dx*dx*dx
}

/// B-spline basis function N_{i,p}(t) using Cox-de Boor recursion.
pub fn bspline_basis(i:usize,p:usize,t:f64,knots:&[f64])->f64{
    if p==0{return if t>=knots[i]&&t<knots[i+1]{1.0}else{0.0};}
    let mut result=0.0;
    let d1=knots[i+p]-knots[i];if d1>1e-30{result+=(t-knots[i])/d1*bspline_basis(i,p-1,t,knots);}
    let d2=knots[i+p+1]-knots[i+1];if d2>1e-30{result+=(knots[i+p+1]-t)/d2*bspline_basis(i+1,p-1,t,knots);}
    result
}

/// Evaluate B-spline curve at parameter t.
pub fn bspline_eval(control_points:&[(f64,f64)],knots:&[f64],degree:usize,t:f64)->(f64,f64){
    let n=control_points.len();let mut x=0.0;let mut y=0.0;
    for i in 0..n{let b=bspline_basis(i,degree,t,knots);x+=b*control_points[i].0;y+=b*control_points[i].1;}(x,y)
}

/// De Boor's algorithm for B-spline evaluation (more numerically stable).
pub fn de_boor(control_points:&[(f64,f64)],knots:&[f64],degree:usize,t:f64)->(f64,f64){
    let n=control_points.len();let p=degree;
    // Find knot span
    let mut k=p;for i in p..n{if t>=knots[i]&&t<knots[i+1]{k=i;break;}}
    let mut d:Vec<(f64,f64)>=(0..=p).map(|j|control_points[k-p+j]).collect();
    for r in 1..=p{for j in (r..=p).rev(){let i=k-p+j;
        let alpha=(t-knots[i])/(knots[i+p+1-r]-knots[i]);
        d[j]=(d[j-1].0*(1.0-alpha)+d[j].0*alpha,d[j-1].1*(1.0-alpha)+d[j].1*alpha);}}
    d[p]
}

/// Generate uniform knot vector for B-spline.
pub fn uniform_knots(n_control:usize,degree:usize)->Vec<f64>{
    let total=n_control+degree+1;
    (0..total).map(|i|i as f64/(total-1)as f64).collect()
}

/// Generate clamped knot vector.
pub fn clamped_knots(n_control:usize,degree:usize)->Vec<f64>{
    let total=n_control+degree+1;let internal=total-2*(degree+1);
    let mut knots=vec![0.0;degree+1];
    for i in 1..=internal{knots.push(i as f64/(internal+1)as f64);}
    knots.extend(vec![1.0;degree+1]);knots
}

/// Cubic Hermite spline interpolation.
pub fn hermite_spline(x:&[f64],y:&[f64],dy:&[f64],t:f64)->f64{
    let n=x.len()-1;let mut i=0;for k in 0..n{if t>=x[k]{i=k;}}
    let h=x[i+1]-x[i];let s=(t-x[i])/h;let s2=s*s;let s3=s2*s;
    let h00=2.0*s3-3.0*s2+1.0;let h10=s3-2.0*s2+s;let h01=-2.0*s3+3.0*s2;let h11=s3-s2;
    h00*y[i]+h10*h*dy[i]+h01*y[i+1]+h11*h*dy[i+1]
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_cubic_spline(){let x=[0.0,1.0,2.0,3.0];let y=[0.0,1.0,0.0,1.0];
        let coeffs=cubic_spline(&x,&y);let v=cubic_spline_eval(&x,&coeffs,0.0);assert!(approx_eq(v,0.0,1e-10));}
    #[test]fn test_spline_interpolates(){let x=[0.0,1.0,2.0,3.0];let y=[1.0,2.0,0.0,3.0];
        let coeffs=cubic_spline(&x,&y);for i in 0..x.len()-1{assert!(approx_eq(cubic_spline_eval(&x,&coeffs,x[i]),y[i],1e-10));}}
    #[test]fn test_bspline_partition(){let knots=[0.0,0.0,0.0,1.0,1.0,1.0];
        // Basis functions should sum to 1
        let t=0.5;let sum:f64=(0..3).map(|i|bspline_basis(i,2,t,&knots)).sum();
        assert!(approx_eq(sum,1.0,1e-10));}
}
