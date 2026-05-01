//! Coordinate system conversions.

use crate::trig;
use crate::powers;

pub fn cartesian_to_polar(x:f64,y:f64)->(f64,f64){(powers::hypot(x,y),trig::atan2(y,x))}
pub fn polar_to_cartesian(r:f64,t:f64)->(f64,f64){(r*trig::cos(t),r*trig::sin(t))}

pub fn cartesian_to_cylindrical(x:f64,y:f64,z:f64)->(f64,f64,f64){let(r,t)=cartesian_to_polar(x,y);(r,t,z)}
pub fn cylindrical_to_cartesian(r:f64,t:f64,z:f64)->(f64,f64,f64){let(x,y)=polar_to_cartesian(r,t);(x,y,z)}

pub fn cartesian_to_spherical(x:f64,y:f64,z:f64)->(f64,f64,f64){
    let r=powers::sqrt(x*x+y*y+z*z);if r==0.0{return(0.0,0.0,0.0);}
    (r,trig::acos(z/r),trig::atan2(y,x))
}
pub fn spherical_to_cartesian(r:f64,theta:f64,phi:f64)->(f64,f64,f64){
    let st=trig::sin(theta);(r*st*trig::cos(phi),r*st*trig::sin(phi),r*trig::cos(theta))
}
pub fn spherical_to_cylindrical(r:f64,theta:f64,phi:f64)->(f64,f64,f64){(r*trig::sin(theta),phi,r*trig::cos(theta))}
pub fn cylindrical_to_spherical(rho:f64,phi:f64,z:f64)->(f64,f64,f64){(powers::sqrt(rho*rho+z*z),trig::atan2(rho,z),phi)}

pub fn geo_to_cartesian(lat_deg:f64,lon_deg:f64)->(f64,f64,f64){
    let(la,lo)=(trig::to_radians(lat_deg),trig::to_radians(lon_deg));
    (trig::cos(la)*trig::cos(lo),trig::cos(la)*trig::sin(lo),trig::sin(la))
}

pub fn haversine_distance(lat1:f64,lon1:f64,lat2:f64,lon2:f64,radius:f64)->f64{
    let(la1,lo1,la2,lo2)=(trig::to_radians(lat1),trig::to_radians(lon1),trig::to_radians(lat2),trig::to_radians(lon2));
    let a=trig::haversin(la2-la1)+trig::cos(la1)*trig::cos(la2)*trig::haversin(lo2-lo1);
    2.0*radius*trig::asin(powers::sqrt(a))
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_polar(){let(r,t)=cartesian_to_polar(3.0,4.0);let(x,y)=polar_to_cartesian(r,t);assert!(approx_eq(x,3.0,1e-10));assert!(approx_eq(y,4.0,1e-10));}
    #[test]fn test_spherical(){let(r,t,p)=cartesian_to_spherical(1.0,2.0,3.0);let(x,y,z)=spherical_to_cartesian(r,t,p);assert!(approx_eq(x,1.0,1e-10));assert!(approx_eq(y,2.0,1e-10));assert!(approx_eq(z,3.0,1e-10));}
    #[test]fn test_haversine(){let d=haversine_distance(40.7128,-74.006,51.5074,-0.1278,6371.0);assert!(d>5500.0&&d<5700.0);}
}
