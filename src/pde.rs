//! PDE solvers using finite differences.

/// Solve 1D heat equation u_t = α * u_xx using explicit finite differences.
/// Returns solution matrix [time_step][spatial_point].
pub fn heat_1d(alpha:f64,length:f64,nx:usize,t_end:f64,nt:usize,initial:&dyn Fn(f64)->f64)->Vec<Vec<f64>>{
    let dx=length/(nx-1)as f64;let dt=t_end/nt as f64;let r=alpha*dt/(dx*dx);
    let mut u:Vec<f64>=(0..nx).map(|i|initial(i as f64*dx)).collect();
    let mut result=vec![u.clone()];
    for _ in 0..nt{let mut u_new=u.clone();
        for i in 1..nx-1{u_new[i]=u[i]+r*(u[i+1]-2.0*u[i]+u[i-1]);}
        u=u_new;result.push(u.clone());}
    result
}

/// Solve 1D wave equation u_tt = c² * u_xx.
pub fn wave_1d(c:f64,length:f64,nx:usize,t_end:f64,nt:usize,initial:&dyn Fn(f64)->f64,initial_vel:&dyn Fn(f64)->f64)->Vec<Vec<f64>>{
    let dx=length/(nx-1)as f64;let dt=t_end/nt as f64;let r2=(c*dt/dx)*(c*dt/dx);
    let mut u_prev:Vec<f64>=(0..nx).map(|i|initial(i as f64*dx)).collect();
    let mut u_curr=u_prev.clone();
    // First step using initial velocity
    for i in 1..nx-1{u_curr[i]=u_prev[i]+dt*initial_vel(i as f64*dx)+0.5*r2*(u_prev[i+1]-2.0*u_prev[i]+u_prev[i-1]);}
    let mut result=vec![u_prev.clone(),u_curr.clone()];
    for _ in 2..=nt{let mut u_next=vec![0.0;nx];
        for i in 1..nx-1{u_next[i]=2.0*u_curr[i]-u_prev[i]+r2*(u_curr[i+1]-2.0*u_curr[i]+u_curr[i-1]);}
        u_prev=u_curr;u_curr=u_next;result.push(u_curr.clone());}
    result
}

/// Solve Laplace equation ∇²u = 0 on 2D grid using Gauss-Seidel iteration.
pub fn laplace_2d(nx:usize,ny:usize,boundary:&dyn Fn(usize,usize,usize,usize)->f64,tol:f64,max_iter:usize)->Vec<Vec<f64>>{
    let mut u=vec![vec![0.0;ny];nx];
    // Set boundary conditions
    for i in 0..nx{u[i][0]=boundary(i,0,nx,ny);u[i][ny-1]=boundary(i,ny-1,nx,ny);}
    for j in 0..ny{u[0][j]=boundary(0,j,nx,ny);u[nx-1][j]=boundary(nx-1,j,nx,ny);}
    for _ in 0..max_iter{let mut max_diff=0.0f64;
        for i in 1..nx-1{for j in 1..ny-1{
            let new_val=0.25*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]);
            let diff=crate::arithmetic::abs(new_val-u[i][j]);if diff>max_diff{max_diff=diff;}
            u[i][j]=new_val;}}
        if max_diff<tol{break;}}
    u
}

/// Solve Poisson equation ∇²u = f on 2D grid.
pub fn poisson_2d(nx:usize,ny:usize,source:&dyn Fn(f64,f64)->f64,dx:f64,dy:f64,tol:f64,max_iter:usize)->Vec<Vec<f64>>{
    let mut u=vec![vec![0.0;ny];nx];
    for _ in 0..max_iter{let mut max_diff=0.0f64;
        for i in 1..nx-1{for j in 1..ny-1{
            let f_val=source(i as f64*dx,j as f64*dy);
            let new_val=0.25*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]-dx*dy*f_val);
            let diff=crate::arithmetic::abs(new_val-u[i][j]);if diff>max_diff{max_diff=diff;}
            u[i][j]=new_val;}}
        if max_diff<tol{break;}}
    u
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_heat(){let r=heat_1d(0.01,1.0,11,0.1,100,&|x|if x>0.4&&x<0.6{1.0}else{0.0});
        assert!(r.len()>1);assert_eq!(r[0].len(),11);}
    #[test]fn test_wave(){let r=wave_1d(1.0,1.0,21,0.5,100,&|x|crate::trig::sin(crate::constants::PI*x),&|_|0.0);
        assert!(r.len()>1);}
    #[test]fn test_laplace(){let u=laplace_2d(11,11,&|i,j,nx,ny|{if j==ny-1{100.0}else{0.0}},1e-6,10000);
        assert!(u[5][5]>0.0&&u[5][5]<100.0);} // interior should be between boundary values
}
