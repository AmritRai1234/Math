//! ODE solvers.

/// Euler method. Returns Vec of (t, y) pairs.
pub fn euler(f:&dyn Fn(f64,f64)->f64, y0:f64, t0:f64, t_end:f64, dt:f64)->Vec<(f64,f64)>{
    let mut result=vec![(t0,y0)];let mut t=t0;let mut y=y0;
    while t<t_end{y+=dt*f(t,y);t+=dt;result.push((t,y));}result
}

/// Midpoint method.
pub fn midpoint(f:&dyn Fn(f64,f64)->f64, y0:f64, t0:f64, t_end:f64, dt:f64)->Vec<(f64,f64)>{
    let mut result=vec![(t0,y0)];let mut t=t0;let mut y=y0;
    while t<t_end{let k1=f(t,y);let k2=f(t+dt/2.0,y+dt/2.0*k1);y+=dt*k2;t+=dt;result.push((t,y));}result
}

/// 4th-order Runge-Kutta.
pub fn rk4(f:&dyn Fn(f64,f64)->f64, y0:f64, t0:f64, t_end:f64, dt:f64)->Vec<(f64,f64)>{
    let mut result=vec![(t0,y0)];let mut t=t0;let mut y=y0;
    while t<t_end{
        let k1=dt*f(t,y);let k2=dt*f(t+dt/2.0,y+k1/2.0);
        let k3=dt*f(t+dt/2.0,y+k2/2.0);let k4=dt*f(t+dt,y+k3);
        y+=(k1+2.0*k2+2.0*k3+k4)/6.0;t+=dt;result.push((t,y));
    }result
}

/// Adaptive RK45 (Dormand-Prince). Returns Vec of (t, y).
pub fn rk45(f:&dyn Fn(f64,f64)->f64, y0:f64, t0:f64, t_end:f64, tol:f64)->Vec<(f64,f64)>{
    let mut result=vec![(t0,y0)];let mut t=t0;let mut y=y0;let mut dt=(t_end-t0)/100.0;
    while t<t_end{
        if t+dt>t_end{dt=t_end-t;}
        let k1=dt*f(t,y);let k2=dt*f(t+dt/4.0,y+k1/4.0);
        let k3=dt*f(t+3.0*dt/8.0,y+3.0*k1/32.0+9.0*k2/32.0);
        let k4=dt*f(t+12.0*dt/13.0,y+1932.0*k1/2197.0-7200.0*k2/2197.0+7296.0*k3/2197.0);
        let k5=dt*f(t+dt,y+439.0*k1/216.0-8.0*k2+3680.0*k3/513.0-845.0*k4/4104.0);
        let k6=dt*f(t+dt/2.0,y-8.0*k1/27.0+2.0*k2-3544.0*k3/2565.0+1859.0*k4/4104.0-11.0*k5/40.0);
        let y4=y+25.0*k1/216.0+1408.0*k3/2565.0+2197.0*k4/4104.0-k5/5.0;
        let y5=y+16.0*k1/135.0+6656.0*k3/12825.0+28561.0*k4/56430.0-9.0*k5/50.0+2.0*k6/55.0;
        let err=crate::arithmetic::abs(y5-y4);
        if err<tol{t+=dt;y=y5;result.push((t,y));}
        let s=if err>1e-30{0.84*(tol/err).powf(0.25)}else{2.0};
        dt*=s.max(0.1).min(4.0);
    }result
}

/// Adams-Bashforth 4-step method.
pub fn adams_bashforth(f:&dyn Fn(f64,f64)->f64, y0:f64, t0:f64, t_end:f64, dt:f64)->Vec<(f64,f64)>{
    let init=rk4(f,y0,t0,(t0+3.0*dt).min(t_end),dt);
    let mut result=init.clone();if init.len()<4{return result;}
    let mut fs:Vec<f64>=init.iter().map(|&(t,y)|f(t,y)).collect();
    let mut t=init.last().unwrap().0;let mut y=init.last().unwrap().1;
    while t<t_end{let n=fs.len();
        y+=dt/24.0*(55.0*fs[n-1]-59.0*fs[n-2]+37.0*fs[n-3]-9.0*fs[n-4]);
        t+=dt;fs.push(f(t,y));result.push((t,y));
    }result
}

/// Velocity Verlet for second-order ODEs: x'' = a(t, x).
pub fn verlet(a:&dyn Fn(f64,f64)->f64, x0:f64, v0:f64, t0:f64, t_end:f64, dt:f64)->Vec<(f64,f64,f64)>{
    let mut result=vec![(t0,x0,v0)];let mut t=t0;let mut x=x0;let mut v=v0;let mut acc=a(t,x);
    while t<t_end{x+=v*dt+0.5*acc*dt*dt;let new_acc=a(t+dt,x);v+=0.5*(acc+new_acc)*dt;acc=new_acc;t+=dt;result.push((t,x,v));}
    result
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_euler(){let r=euler(&|_,y|y,1.0,0.0,1.0,0.001);let last=r.last().unwrap().1;assert!(approx_eq(last,crate::constants::E,0.01));}
    #[test]fn test_rk4(){let r=rk4(&|_,y|y,1.0,0.0,1.0,0.01);let last=r.last().unwrap().1;assert!(approx_eq(last,crate::constants::E,1e-6));}
}
