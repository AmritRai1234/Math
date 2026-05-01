//! Dynamical systems: Lorenz, Rössler, logistic map, Lyapunov exponents.

use crate::powers;use crate::exponential;

/// Lorenz attractor: dx/dt = σ(y-x), dy/dt = x(ρ-z)-y, dz/dt = xy-βz.
pub fn lorenz(x0:f64,y0:f64,z0:f64,sigma:f64,rho:f64,beta:f64,dt:f64,steps:usize)->Vec<(f64,f64,f64)>{
    let(mut x,mut y,mut z)=(x0,y0,z0);let mut result=vec![(x,y,z)];
    for _ in 0..steps{let dx=sigma*(y-x)*dt;let dy=(x*(rho-z)-y)*dt;let dz=(x*y-beta*z)*dt;
        x+=dx;y+=dy;z+=dz;result.push((x,y,z));}result
}

/// Rössler attractor: dx=-y-z, dy=x+ay, dz=b+z(x-c).
pub fn rossler(x0:f64,y0:f64,z0:f64,a:f64,b:f64,c:f64,dt:f64,steps:usize)->Vec<(f64,f64,f64)>{
    let(mut x,mut y,mut z)=(x0,y0,z0);let mut result=vec![(x,y,z)];
    for _ in 0..steps{let dx=(-y-z)*dt;let dy=(x+a*y)*dt;let dz=(b+z*(x-c))*dt;
        x+=dx;y+=dy;z+=dz;result.push((x,y,z));}result
}

/// Logistic map: x_{n+1} = r * x_n * (1 - x_n).
pub fn logistic_map(r:f64,x0:f64,steps:usize)->Vec<f64>{
    let mut x=x0;let mut result=vec![x];
    for _ in 0..steps{x=r*x*(1.0-x);result.push(x);}result
}

/// Hénon map: x_{n+1} = 1 - a*x_n² + y_n, y_{n+1} = b*x_n.
pub fn henon_map(a:f64,b:f64,x0:f64,y0:f64,steps:usize)->Vec<(f64,f64)>{
    let(mut x,mut y)=(x0,y0);let mut result=vec![(x,y)];
    for _ in 0..steps{let xn=1.0-a*x*x+y;let yn=b*x;x=xn;y=yn;result.push((x,y));}result
}

/// Lyapunov exponent of logistic map.
pub fn lyapunov_logistic(r:f64,x0:f64,transient:usize,iterations:usize)->f64{
    let mut x=x0;for _ in 0..transient{x=r*x*(1.0-x);}
    let mut sum=0.0;for _ in 0..iterations{x=r*x*(1.0-x);
        let deriv=crate::arithmetic::abs(r*(1.0-2.0*x));
        if deriv>0.0{sum+=exponential::ln(deriv);}}
    sum/iterations as f64
}

/// Bifurcation diagram data: for each r, returns the attractor values.
pub fn bifurcation_diagram(r_start:f64,r_end:f64,n_r:usize,transient:usize,collect:usize)->Vec<(f64,f64)>{
    let mut data=Vec::new();
    for i in 0..n_r{let r=r_start+(r_end-r_start)*i as f64/n_r as f64;
        let mut x=0.5;for _ in 0..transient{x=r*x*(1.0-x);}
        for _ in 0..collect{x=r*x*(1.0-x);data.push((r,x));}}data
}

/// Feigenbaum constant δ ≈ 4.669... (ratio of bifurcation widths).
pub fn feigenbaum_delta()->f64{4.669201609102990}

/// Tent map: x_{n+1} = μ*min(x, 1-x).
pub fn tent_map(mu:f64,x0:f64,steps:usize)->Vec<f64>{
    let mut x=x0;let mut result=vec![x];
    for _ in 0..steps{x=mu*x.min(1.0-x);result.push(x);}result
}

/// Baker's map.
pub fn bakers_map(x:f64,y:f64)->(f64,f64){
    if x<0.5{(2.0*x,y/2.0)}else{(2.0*x-1.0,(y+1.0)/2.0)}
}

/// Arnold's cat map on [0,1)×[0,1).
pub fn cat_map(x:f64,y:f64)->(f64,f64){((2.0*x+y)%1.0,(x+y)%1.0)}

/// Standard map (Chirikov).
pub fn standard_map(theta:f64,p:f64,k:f64)->(f64,f64){
    let pn=p+k*crate::trig::sin(theta)/(2.0*crate::constants::PI);
    let tn=theta+pn;(tn%1.0,(pn%1.0+1.0)%1.0)
}

/// Detect period of iterated map.
pub fn detect_period(values:&[f64],tol:f64)->Option<usize>{
    let n=values.len();if n<10{return None;}
    for period in 1..n/2{let last=values[n-1];
        if(values[n-1-period]-last).abs()<tol{
            let mut valid=true;for k in 1..5.min(n/period){
                if(values[n-1-k*period]-last).abs()>tol{valid=false;break;}}
            if valid{return Some(period);}}}None
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_lorenz_starts(){let traj=lorenz(1.0,1.0,1.0,10.0,28.0,8.0/3.0,0.001,100);
        assert_eq!(traj.len(),101);}
    #[test]fn test_logistic_fixed(){let seq=logistic_map(2.0,0.5,100);// r=2: fixed point at 0.5
        assert!(approx_eq(*seq.last().unwrap(),0.5,1e-6));}
    #[test]fn test_lyapunov_chaos(){let le=lyapunov_logistic(4.0,0.1,1000,10000);assert!(le>0.0);} // chaotic
    #[test]fn test_feigenbaum(){assert!(approx_eq(feigenbaum_delta(),4.669,0.001));}
}
