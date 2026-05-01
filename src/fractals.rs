//! Fractals: Mandelbrot, Julia, IFS, fractal dimension.

use crate::complex::Complex;

/// Mandelbrot iteration count for point c. Returns iterations before escape (max if in set).
pub fn mandelbrot(cr:f64,ci:f64,max_iter:u32)->u32{
    let mut zr=0.0;let mut zi=0.0;
    for i in 0..max_iter{let zr2=zr*zr;let zi2=zi*zi;
        if zr2+zi2>4.0{return i;}zi=2.0*zr*zi+ci;zr=zr2-zi2+cr;}max_iter
}

/// Julia set iteration count for point z with constant c.
pub fn julia(zr:f64,zi:f64,cr:f64,ci:f64,max_iter:u32)->u32{
    let(mut zr,mut zi)=(zr,zi);
    for i in 0..max_iter{let zr2=zr*zr;let zi2=zi*zi;
        if zr2+zi2>4.0{return i;}let new_zr=zr2-zi2+cr;zi=2.0*zr*zi+ci;zr=new_zr;}max_iter
}

/// Generate Mandelbrot set as a grid of iteration counts.
pub fn mandelbrot_grid(x_min:f64,x_max:f64,y_min:f64,y_max:f64,width:usize,height:usize,max_iter:u32)->Vec<Vec<u32>>{
    let dx=(x_max-x_min)/width as f64;let dy=(y_max-y_min)/height as f64;
    (0..height).map(|j|{let ci=y_min+j as f64*dy;
        (0..width).map(|i|mandelbrot(x_min+i as f64*dx,ci,max_iter)).collect()}).collect()
}

/// Burning Ship fractal.
pub fn burning_ship(cr:f64,ci:f64,max_iter:u32)->u32{
    let(mut zr,mut zi)=(0.0,0.0);
    for i in 0..max_iter{let zr2=zr*zr;let zi2=zi*zi;
        if zr2+zi2>4.0{return i;}let new_zr=zr2-zi2+cr;
        zi=crate::arithmetic::abs(2.0*zr*zi)+ci;zr=new_zr;}max_iter
}

/// Sierpinski triangle: check if point (x,y) is in the triangle at given depth.
pub fn sierpinski(x:u32,y:u32)->bool{(x&y)==0}

/// Cantor set: check if a number (in ternary) avoids digit 1.
pub fn in_cantor_set(mut n:u64,base:u64)->bool{
    while n>0{if n%base==1{return false;}n/=base;}true
}

/// Koch snowflake perimeter after n iterations (starting side length s).
pub fn koch_perimeter(s:f64,n:u32)->f64{s*3.0*crate::powers::pow(4.0/3.0,n as f64)}

/// Menger sponge volume after n iterations (starting side s).
pub fn menger_volume(s:f64,n:u32)->f64{crate::powers::pow(s,3.0)*crate::powers::pow(20.0/27.0,n as f64)}

/// Fractal dimension estimate using box-counting.
pub fn box_counting_dimension(sizes:&[f64],counts:&[f64])->f64{
    let n=sizes.len();if n<2{return 0.0;}
    let log_sizes:Vec<f64>=sizes.iter().map(|&s|crate::exponential::ln(1.0/s)).collect();
    let log_counts:Vec<f64>=counts.iter().map(|&c|crate::exponential::ln(c)).collect();
    let(_,slope)=crate::regression::linear_regression(&log_sizes,&log_counts);slope
}

/// Logistic map: x_{n+1} = r * x_n * (1 - x_n). Returns orbit.
pub fn logistic_map(r:f64,x0:f64,n:usize)->Vec<f64>{
    let mut orbit=vec![x0];let mut x=x0;
    for _ in 0..n{x=r*x*(1.0-x);orbit.push(x);}orbit
}

/// Henon map: returns orbit of (x, y) pairs.
pub fn henon_map(a:f64,b:f64,x0:f64,y0:f64,n:usize)->Vec<(f64,f64)>{
    let mut orbit=vec![(x0,y0)];let(mut x,mut y)=(x0,y0);
    for _ in 0..n{let new_x=1.0-a*x*x+y;let new_y=b*x;x=new_x;y=new_y;orbit.push((x,y));}orbit
}

/// Lyapunov exponent for the logistic map.
pub fn lyapunov_exponent(r:f64,x0:f64,n:usize,skip:usize)->f64{
    let mut x=x0;for _ in 0..skip{x=r*x*(1.0-x);}
    let mut sum=0.0;for _ in 0..n{let deriv=crate::arithmetic::abs(r-2.0*r*x);
        if deriv>1e-30{sum+=crate::exponential::ln(deriv);}x=r*x*(1.0-x);}sum/n as f64
}

/// Iterated function system (IFS): apply random affine transforms.
pub fn ifs_iterate(transforms:&[(f64,f64,f64,f64,f64,f64)],weights:&[f64],n:usize,seed:u64)->Vec<(f64,f64)>{
    let mut rng=seed;let(mut x,mut y)=(0.0,0.0);let mut points=Vec::with_capacity(n);
    let total:f64=weights.iter().sum();
    for _ in 0..n{rng=crate::noise::xorshift(rng);let r=(rng%10000)as f64/10000.0*total;
        let mut cum=0.0;let mut idx=0;for(i,&w) in weights.iter().enumerate(){cum+=w;if r<=cum{idx=i;break;}}
        let(a,b,c,d,e,f)=transforms[idx];let nx=a*x+b*y+e;let ny=c*x+d*y+f;x=nx;y=ny;
        points.push((x,y));}points
}

/// Barnsley fern IFS transforms.
pub fn barnsley_fern(n:usize,seed:u64)->Vec<(f64,f64)>{
    let transforms=[(0.0,0.0,0.0,0.16,0.0,0.0),(0.85,0.04,-0.04,0.85,0.0,1.6),
        (0.2,-0.26,0.23,0.22,0.0,1.6),(-0.15,0.28,0.26,0.24,0.0,0.44)];
    let weights=[0.01,0.85,0.07,0.07];
    ifs_iterate(&transforms,&weights,n,seed)
}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_mandelbrot_in(){assert_eq!(mandelbrot(0.0,0.0,100),100);} // origin is in set
    #[test]fn test_mandelbrot_out(){assert!(mandelbrot(2.0,2.0,100)<100);} // far point escapes
    #[test]fn test_julia(){assert!(julia(0.0,0.0,0.0,0.0,100)==100);} // c=0 → z stays at 0
    #[test]fn test_sierpinski(){assert!(sierpinski(0,0));assert!(!sierpinski(3,5));}
    #[test]fn test_logistic(){let orbit=logistic_map(3.2,0.5,100);assert_eq!(orbit.len(),101);}
    #[test]fn test_koch(){let p=koch_perimeter(1.0,0);assert!((p-3.0).abs()<1e-10);}
}
