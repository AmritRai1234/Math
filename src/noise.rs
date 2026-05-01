//! Noise generation: Perlin, simplex, fBm, RNG.

use crate::arithmetic;
use crate::interpolation;

// ─── Simple RNG ───

/// Linear congruential generator.
pub fn lcg(seed:u64)->u64{seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407)}

/// Xorshift64 RNG.
pub fn xorshift(mut seed:u64)->u64{seed^=seed<<13;seed^=seed>>7;seed^=seed<<17;seed}

/// White noise: returns f64 in [0, 1) from seed.
pub fn white_noise(seed:u64)->f64{(xorshift(seed)>>11)as f64/9007199254740992.0}

// ─── Perlin Noise ───

fn fade(t:f64)->f64{t*t*t*(t*(t*6.0-15.0)+10.0)}

fn hash2d(x:i32,y:i32)->u64{
    let mut h=(x as u64).wrapping_mul(374761393).wrapping_add((y as u64).wrapping_mul(668265263));
    h=(h^(h>>13)).wrapping_mul(1274126177);h^(h>>16)
}

fn grad2d(hash:u64,x:f64,y:f64)->f64{
    match hash&3{0=>x+y,1=>-x+y,2=>x-y,_=>-x-y}
}

/// 2D Perlin noise. Returns value in approximately [-1, 1].
pub fn perlin_noise_2d(x:f64,y:f64)->f64{
    let xi=arithmetic::floor(x)as i32;let yi=arithmetic::floor(y)as i32;
    let xf=x-xi as f64;let yf=y-yi as f64;
    let u=fade(xf);let v=fade(yf);
    let n00=grad2d(hash2d(xi,yi),xf,yf);
    let n10=grad2d(hash2d(xi+1,yi),xf-1.0,yf);
    let n01=grad2d(hash2d(xi,yi+1),xf,yf-1.0);
    let n11=grad2d(hash2d(xi+1,yi+1),xf-1.0,yf-1.0);
    let nx0=interpolation::lerp(n00,n10,u);
    let nx1=interpolation::lerp(n01,n11,u);
    interpolation::lerp(nx0,nx1,v)
}

fn hash3d(x:i32,y:i32,z:i32)->u64{
    let mut h=(x as u64).wrapping_mul(374761393).wrapping_add((y as u64).wrapping_mul(668265263)).wrapping_add((z as u64).wrapping_mul(1440670441));
    h=(h^(h>>13)).wrapping_mul(1274126177);h^(h>>16)
}

fn grad3d(hash:u64,x:f64,y:f64,z:f64)->f64{
    match hash%12{0=>x+y,1=>-x+y,2=>x-y,3=>-x-y,4=>x+z,5=>-x+z,6=>x-z,7=>-x-z,8=>y+z,9=>-y+z,10=>y-z,_=>-y-z}
}

/// 3D Perlin noise.
pub fn perlin_noise_3d(x:f64,y:f64,z:f64)->f64{
    let xi=arithmetic::floor(x)as i32;let yi=arithmetic::floor(y)as i32;let zi=arithmetic::floor(z)as i32;
    let xf=x-xi as f64;let yf=y-yi as f64;let zf=z-zi as f64;
    let u=fade(xf);let v=fade(yf);let w=fade(zf);
    let n000=grad3d(hash3d(xi,yi,zi),xf,yf,zf);let n100=grad3d(hash3d(xi+1,yi,zi),xf-1.0,yf,zf);
    let n010=grad3d(hash3d(xi,yi+1,zi),xf,yf-1.0,zf);let n110=grad3d(hash3d(xi+1,yi+1,zi),xf-1.0,yf-1.0,zf);
    let n001=grad3d(hash3d(xi,yi,zi+1),xf,yf,zf-1.0);let n101=grad3d(hash3d(xi+1,yi,zi+1),xf-1.0,yf,zf-1.0);
    let n011=grad3d(hash3d(xi,yi+1,zi+1),xf,yf-1.0,zf-1.0);let n111=grad3d(hash3d(xi+1,yi+1,zi+1),xf-1.0,yf-1.0,zf-1.0);
    let nx00=interpolation::lerp(n000,n100,u);let nx10=interpolation::lerp(n010,n110,u);
    let nx01=interpolation::lerp(n001,n101,u);let nx11=interpolation::lerp(n011,n111,u);
    let nxy0=interpolation::lerp(nx00,nx10,v);let nxy1=interpolation::lerp(nx01,nx11,v);
    interpolation::lerp(nxy0,nxy1,w)
}

/// 2D Simplex noise (simplified).
pub fn simplex_noise_2d(x:f64,y:f64)->f64{
    // Skew to simplex space
    let f2=0.5*(crate::powers::sqrt(3.0)-1.0);let g2=(3.0-crate::powers::sqrt(3.0))/6.0;
    let s=(x+y)*f2;let i=arithmetic::floor(x+s)as i32;let j=arithmetic::floor(y+s)as i32;
    let t=(i+j)as f64*g2;let x0=x-(i as f64-t);let y0=y-(j as f64-t);
    let(i1,j1)=if x0>y0{(1,0)}else{(0,1)};
    let x1=x0-i1 as f64+g2;let y1=y0-j1 as f64+g2;let x2=x0-1.0+2.0*g2;let y2=y0-1.0+2.0*g2;
    let mut n0=0.0;let t0=0.5-x0*x0-y0*y0;if t0>0.0{let t0=t0*t0;n0=t0*t0*grad2d(hash2d(i,j),x0,y0);}
    let mut n1=0.0;let t1=0.5-x1*x1-y1*y1;if t1>0.0{let t1=t1*t1;n1=t1*t1*grad2d(hash2d(i+i1,j+j1),x1,y1);}
    let mut n2=0.0;let t2=0.5-x2*x2-y2*y2;if t2>0.0{let t2=t2*t2;n2=t2*t2*grad2d(hash2d(i+1,j+1),x2,y2);}
    70.0*(n0+n1+n2)
}

/// Fractal Brownian Motion.
pub fn fbm(x:f64,y:f64,octaves:u32,lacunarity:f64,gain:f64)->f64{
    let mut sum=0.0;let mut amp=1.0;let mut freq=1.0;
    for _ in 0..octaves{sum+=amp*perlin_noise_2d(x*freq,y*freq);freq*=lacunarity;amp*=gain;}sum
}

#[cfg(test)]
mod tests{
    use super::*;
    #[test]fn test_perlin_range(){let v=perlin_noise_2d(0.5,0.5);assert!(v>=-2.0&&v<=2.0);}
    #[test]fn test_perlin_deterministic(){assert_eq!(perlin_noise_2d(1.0,1.0),perlin_noise_2d(1.0,1.0));}
    #[test]fn test_lcg(){let a=lcg(42);let b=lcg(42);assert_eq!(a,b);assert_ne!(a,42);}
}
