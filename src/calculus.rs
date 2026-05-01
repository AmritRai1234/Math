//! Numerical calculus: derivatives, integrals, gradient, Jacobian, Hessian.

/// Numerical derivative using central difference.
pub fn derivative(f:&dyn Fn(f64)->f64, x:f64, h:f64)->f64{(f(x+h)-f(x-h))/(2.0*h)}

/// Second derivative.
pub fn second_derivative(f:&dyn Fn(f64)->f64, x:f64, h:f64)->f64{(f(x+h)-2.0*f(x)+f(x-h))/(h*h)}

/// Partial derivative of f(point) with respect to variable `var`.
pub fn partial_deriv(f:&dyn Fn(&[f64])->f64, point:&[f64], var:usize, h:f64)->f64{
    let mut p1=point.to_vec();let mut p2=point.to_vec();
    p1[var]+=h;p2[var]-=h;(f(&p1)-f(&p2))/(2.0*h)
}

/// Gradient vector.
pub fn gradient(f:&dyn Fn(&[f64])->f64, point:&[f64], h:f64)->Vec<f64>{
    (0..point.len()).map(|i|partial_deriv(f,point,i,h)).collect()
}

/// Jacobian matrix for a system of functions.
pub fn jacobian(fs:&[&dyn Fn(&[f64])->f64], point:&[f64], h:f64)->Vec<Vec<f64>>{
    fs.iter().map(|fi|(0..point.len()).map(|j|partial_deriv(fi,point,j,h)).collect()).collect()
}

/// Hessian matrix (second partial derivatives).
pub fn hessian(f:&dyn Fn(&[f64])->f64, point:&[f64], h:f64)->Vec<Vec<f64>>{
    let n=point.len();let mut result=vec![vec![0.0;n];n];
    for i in 0..n{for j in 0..n{
        let mut pp=point.to_vec();let mut pm=point.to_vec();let mut mp=point.to_vec();let mut mm=point.to_vec();
        pp[i]+=h;pp[j]+=h;pm[i]+=h;pm[j]-=h;mp[i]-=h;mp[j]+=h;mm[i]-=h;mm[j]-=h;
        result[i][j]=(f(&pp)-f(&pm)-f(&mp)+f(&mm))/(4.0*h*h);
    }}
    result
}

/// Trapezoidal rule integration.
pub fn integrate_trapezoid(f:&dyn Fn(f64)->f64, a:f64, b:f64, n:usize)->f64{
    let h=(b-a)/n as f64;let mut sum=0.5*(f(a)+f(b));
    for i in 1..n{sum+=f(a+i as f64*h);}sum*h
}

/// Simpson's rule integration.
pub fn integrate_simpson(f:&dyn Fn(f64)->f64, a:f64, b:f64, n:usize)->f64{
    let n=if n%2!=0{n+1}else{n};let h=(b-a)/n as f64;
    let mut sum=f(a)+f(b);
    for i in 1..n{let x=a+i as f64*h;sum+=if i%2==0{2.0*f(x)}else{4.0*f(x)};}
    sum*h/3.0
}

/// Romberg integration.
pub fn integrate_romberg(f:&dyn Fn(f64)->f64, a:f64, b:f64, n:usize)->f64{
    let mut r=vec![vec![0.0;n];n];r[0][0]=(b-a)/2.0*(f(a)+f(b));
    for i in 1..n{let h=(b-a)/((1<<i)as f64);let mut sum=0.0;
        for k in 1..=(1<<(i-1)){sum+=f(a+(2*k-1)as f64*h);}
        r[i][0]=0.5*r[i-1][0]+h*sum;
        for j in 1..=i{let p=crate::powers::powi(4.0,j as u32);r[i][j]=(p*r[i][j-1]-r[i-1][j-1])/(p-1.0);}
    }
    r[n-1][n-1]
}

/// Gaussian quadrature (5-point Gauss-Legendre).
pub fn integrate_gauss(f:&dyn Fn(f64)->f64, a:f64, b:f64, _n:usize)->f64{
    let nodes=[0.0,-0.5384693101056831,0.5384693101056831,-0.9061798459386640,0.9061798459386640];
    let weights=[0.5688888888888889,0.4786286704993665,0.4786286704993665,0.2369268850561891,0.2369268850561891];
    let mid=(a+b)/2.0;let half=(b-a)/2.0;
    let mut sum=0.0;for i in 0..5{sum+=weights[i]*f(mid+half*nodes[i]);}sum*half
}

/// Monte Carlo integration.
pub fn integrate_monte_carlo(f:&dyn Fn(f64)->f64, a:f64, b:f64, n:usize)->f64{
    let range=b-a;let mut sum=0.0;let mut rng=12345u64;
    for _ in 0..n{rng=rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        let t=(rng>>33)as f64/2147483648.0;sum+=f(a+t*range);}
    sum*range/n as f64
}

/// Divergence of a 3D vector field F at a point.
pub fn divergence(f:&[&dyn Fn(&[f64])->f64;3], point:&[f64], h:f64)->f64{
    partial_deriv(f[0],point,0,h)+partial_deriv(f[1],point,1,h)+partial_deriv(f[2],point,2,h)
}

/// Curl of a 3D vector field.
pub fn curl(f:&[&dyn Fn(&[f64])->f64;3], point:&[f64], h:f64)->[f64;3]{[
    partial_deriv(f[2],point,1,h)-partial_deriv(f[1],point,2,h),
    partial_deriv(f[0],point,2,h)-partial_deriv(f[2],point,0,h),
    partial_deriv(f[1],point,0,h)-partial_deriv(f[0],point,1,h),
]}

/// Laplacian of a scalar field.
pub fn laplacian(f:&dyn Fn(&[f64])->f64, point:&[f64], h:f64)->f64{
    let n=point.len();let mut sum=0.0;
    for i in 0..n{let mut pp=point.to_vec();let mut pm=point.to_vec();pp[i]+=h;pm[i]-=h;
        sum+=(f(&pp)-2.0*f(point)+f(&pm))/(h*h);}
    sum
}

#[cfg(test)]
mod tests{
    use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_derivative(){let f=|x:f64|x*x;assert!(approx_eq(derivative(&f,3.0,1e-6),6.0,1e-4));}
    #[test]fn test_simpson(){let f=|x:f64|x*x;assert!(approx_eq(integrate_simpson(&f,0.0,1.0,100),1.0/3.0,1e-8));}
    #[test]fn test_gradient(){let f=|p:&[f64]|p[0]*p[0]+p[1]*p[1];let g=gradient(&f,&[3.0,4.0],1e-6);assert!(approx_eq(g[0],6.0,1e-4));assert!(approx_eq(g[1],8.0,1e-4));}
}
