//! Coordinate system conversions: polar, cylindrical, spherical, toroidal.

use crate::trig;use crate::powers;use crate::constants::PI;

// ─── 2D ───

/// Cartesian (x,y) → Polar (r, θ).
pub fn cart_to_polar(x:f64,y:f64)->(f64,f64){(powers::sqrt(x*x+y*y),trig::atan2(y,x))}
/// Polar → Cartesian.
pub fn polar_to_cart(r:f64,theta:f64)->(f64,f64){(r*trig::cos(theta),r*trig::sin(theta))}

// ─── 3D ───

/// Cartesian → Cylindrical (r, θ, z).
pub fn cart_to_cylindrical(x:f64,y:f64,z:f64)->(f64,f64,f64){let(r,t)=cart_to_polar(x,y);(r,t,z)}
/// Cylindrical → Cartesian.
pub fn cylindrical_to_cart(r:f64,theta:f64,z:f64)->(f64,f64,f64){let(x,y)=polar_to_cart(r,theta);(x,y,z)}

/// Cartesian → Spherical (r, θ, φ) where θ=polar angle, φ=azimuthal.
pub fn cart_to_spherical(x:f64,y:f64,z:f64)->(f64,f64,f64){
    let r=powers::sqrt(x*x+y*y+z*z);
    let theta=if r>0.0{trig::acos(z/r)}else{0.0};
    let phi=trig::atan2(y,x);(r,theta,phi)
}
/// Spherical → Cartesian.
pub fn spherical_to_cart(r:f64,theta:f64,phi:f64)->(f64,f64,f64){
    (r*trig::sin(theta)*trig::cos(phi),r*trig::sin(theta)*trig::sin(phi),r*trig::cos(theta))
}

/// Cylindrical → Spherical.
pub fn cylindrical_to_spherical(r:f64,theta:f64,z:f64)->(f64,f64,f64){
    let(x,y,z2)=cylindrical_to_cart(r,theta,z);cart_to_spherical(x,y,z2)
}
/// Spherical → Cylindrical.
pub fn spherical_to_cylindrical(r:f64,theta:f64,phi:f64)->(f64,f64,f64){
    let(x,y,z)=spherical_to_cart(r,theta,phi);cart_to_cylindrical(x,y,z)
}

// ─── Toroidal ───

/// Toroidal → Cartesian. (R=major radius, r=minor, u=major angle, v=minor angle).
pub fn toroidal_to_cart(big_r:f64,little_r:f64,u:f64,v:f64)->(f64,f64,f64){
    let x=(big_r+little_r*trig::cos(v))*trig::cos(u);
    let y=(big_r+little_r*trig::cos(v))*trig::sin(u);
    let z=little_r*trig::sin(v);(x,y,z)
}

// ─── Oblate spheroidal ───

/// Oblate spheroidal → Cartesian.
pub fn oblate_spheroidal_to_cart(a:f64,mu:f64,nu:f64,phi:f64)->(f64,f64,f64){
    let x=a*crate::hyperbolic::cosh(mu)*trig::cos(nu)*trig::cos(phi);
    let y=a*crate::hyperbolic::cosh(mu)*trig::cos(nu)*trig::sin(phi);
    let z=a*crate::hyperbolic::sinh(mu)*trig::sin(nu);(x,y,z)
}

// ─── Parabolic ───

/// Parabolic cylindrical → Cartesian.
pub fn parabolic_to_cart(u:f64,v:f64,z:f64)->(f64,f64,f64){(0.5*(u*u-v*v),u*v,z)}

// ─── Geographic ───

/// Geodetic (lat, lon, alt) → ECEF (Earth-Centered Earth-Fixed).
pub fn geodetic_to_ecef(lat:f64,lon:f64,alt:f64)->(f64,f64,f64){
    let a=6378137.0;let f=1.0/298.257223563;let e2=2.0*f-f*f;
    let n=a/powers::sqrt(1.0-e2*trig::sin(lat)*trig::sin(lat));
    ((n+alt)*trig::cos(lat)*trig::cos(lon),(n+alt)*trig::cos(lat)*trig::sin(lon),(n*(1.0-e2)+alt)*trig::sin(lat))
}

/// ECEF → Geodetic (iterative).
pub fn ecef_to_geodetic(x:f64,y:f64,z:f64)->(f64,f64,f64){
    let a=6378137.0;let f=1.0/298.257223563;let e2=2.0*f-f*f;
    let lon=trig::atan2(y,x);let p=powers::sqrt(x*x+y*y);
    let mut lat=trig::atan2(z,p*(1.0-e2));
    for _ in 0..10{let n=a/powers::sqrt(1.0-e2*trig::sin(lat)*trig::sin(lat));
        lat=trig::atan2(z+e2*n*trig::sin(lat),p);}
    let n=a/powers::sqrt(1.0-e2*trig::sin(lat)*trig::sin(lat));
    let alt=p/trig::cos(lat)-n;(lat,lon,alt)
}

/// Rotate point (x,y) by angle θ.
pub fn rotate_2d(x:f64,y:f64,theta:f64)->(f64,f64){
    (x*trig::cos(theta)-y*trig::sin(theta),x*trig::sin(theta)+y*trig::cos(theta))
}

/// Rotate point (x,y,z) around Z axis.
pub fn rotate_z(x:f64,y:f64,z:f64,theta:f64)->(f64,f64,f64){let(rx,ry)=rotate_2d(x,y,theta);(rx,ry,z)}
/// Rotate around X axis.
pub fn rotate_x(x:f64,y:f64,z:f64,theta:f64)->(f64,f64,f64){let(ry,rz)=rotate_2d(y,z,theta);(x,ry,rz)}
/// Rotate around Y axis.
pub fn rotate_y(x:f64,y:f64,z:f64,theta:f64)->(f64,f64,f64){let(rz,rx)=rotate_2d(z,x,theta);(rx,y,rz)}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_polar_roundtrip(){let(r,t)=cart_to_polar(3.0,4.0);let(x,y)=polar_to_cart(r,t);
        assert!(approx_eq(x,3.0,1e-10));assert!(approx_eq(y,4.0,1e-10));}
    #[test]fn test_spherical_roundtrip(){let(r,t,p)=cart_to_spherical(1.0,2.0,3.0);let(x,y,z)=spherical_to_cart(r,t,p);
        assert!(approx_eq(x,1.0,1e-10));assert!(approx_eq(y,2.0,1e-10));assert!(approx_eq(z,3.0,1e-10));}
    #[test]fn test_cylindrical_roundtrip(){let(r,t,z)=cart_to_cylindrical(3.0,4.0,5.0);let(x,y,z2)=cylindrical_to_cart(r,t,z);
        assert!(approx_eq(x,3.0,1e-10));assert!(approx_eq(y,4.0,1e-10));assert!(approx_eq(z2,5.0,1e-10));}
    #[test]fn test_rotate(){let(x,y)=rotate_2d(1.0,0.0,PI/2.0);assert!(approx_eq(x,0.0,1e-10));assert!(approx_eq(y,1.0,1e-10));}
}
