//! Quaternion algebra: rotations, SLERP, interpolation.

use crate::powers;use crate::trig;use crate::constants::PI;

#[derive(Clone,Copy,Debug,PartialEq)]
pub struct Quaternion{pub w:f64,pub x:f64,pub y:f64,pub z:f64}

impl Quaternion{
    pub fn new(w:f64,x:f64,y:f64,z:f64)->Self{Self{w,x,y,z}}
    pub fn identity()->Self{Self{w:1.0,x:0.0,y:0.0,z:0.0}}
    pub fn from_axis_angle(axis:&[f64;3],angle:f64)->Self{
        let half=angle/2.0;let s=trig::sin(half);let n=powers::sqrt(axis[0]*axis[0]+axis[1]*axis[1]+axis[2]*axis[2]);
        if n<1e-30{return Self::identity();}
        Self{w:trig::cos(half),x:axis[0]/n*s,y:axis[1]/n*s,z:axis[2]/n*s}
    }
    pub fn from_euler(roll:f64,pitch:f64,yaw:f64)->Self{
        let(cr,sr)=(trig::cos(roll/2.0),trig::sin(roll/2.0));
        let(cp,sp)=(trig::cos(pitch/2.0),trig::sin(pitch/2.0));
        let(cy,sy)=(trig::cos(yaw/2.0),trig::sin(yaw/2.0));
        Self{w:cr*cp*cy+sr*sp*sy,x:sr*cp*cy-cr*sp*sy,y:cr*sp*cy+sr*cp*sy,z:cr*cp*sy-sr*sp*cy}
    }
    pub fn norm(self)->f64{powers::sqrt(self.w*self.w+self.x*self.x+self.y*self.y+self.z*self.z)}
    pub fn normalize(self)->Self{let n=self.norm();Self{w:self.w/n,x:self.x/n,y:self.y/n,z:self.z/n}}
    pub fn conjugate(self)->Self{Self{w:self.w,x:-self.x,y:-self.y,z:-self.z}}
    pub fn inverse(self)->Self{let n2=self.w*self.w+self.x*self.x+self.y*self.y+self.z*self.z;
        Self{w:self.w/n2,x:-self.x/n2,y:-self.y/n2,z:-self.z/n2}}
    pub fn mul(self,r:Self)->Self{Self{
        w:self.w*r.w-self.x*r.x-self.y*r.y-self.z*r.z,
        x:self.w*r.x+self.x*r.w+self.y*r.z-self.z*r.y,
        y:self.w*r.y-self.x*r.z+self.y*r.w+self.z*r.x,
        z:self.w*r.z+self.x*r.y-self.y*r.x+self.z*r.w}}
    pub fn add(self,r:Self)->Self{Self{w:self.w+r.w,x:self.x+r.x,y:self.y+r.y,z:self.z+r.z}}
    pub fn scale(self,s:f64)->Self{Self{w:self.w*s,x:self.x*s,y:self.y*s,z:self.z*s}}
    pub fn dot(self,r:Self)->f64{self.w*r.w+self.x*r.x+self.y*r.y+self.z*r.z}
    /// Rotate a 3D vector.
    pub fn rotate_vector(self,v:&[f64;3])->[f64;3]{
        let q=self.normalize();let p=Quaternion::new(0.0,v[0],v[1],v[2]);
        let r=q.mul(p).mul(q.conjugate());[r.x,r.y,r.z]
    }
    /// To rotation matrix (3×3).
    pub fn to_matrix(self)->[[f64;3];3]{let q=self.normalize();
        let(w,x,y,z)=(q.w,q.x,q.y,q.z);
        [[1.0-2.0*(y*y+z*z),2.0*(x*y-w*z),2.0*(x*z+w*y)],
         [2.0*(x*y+w*z),1.0-2.0*(x*x+z*z),2.0*(y*z-w*x)],
         [2.0*(x*z-w*y),2.0*(y*z+w*x),1.0-2.0*(x*x+y*y)]]}
    /// To axis-angle.
    pub fn to_axis_angle(self)->([f64;3],f64){let q=self.normalize();
        let angle=2.0*trig::acos(q.w.min(1.0).max(-1.0));let s=powers::sqrt(1.0-q.w*q.w);
        if s<1e-10{([1.0,0.0,0.0],angle)}else{([q.x/s,q.y/s,q.z/s],angle)}}
    /// To Euler angles (roll, pitch, yaw).
    pub fn to_euler(self)->(f64,f64,f64){let q=self.normalize();
        let sinp=2.0*(q.w*q.y-q.z*q.x);
        let pitch=if sinp.abs()>=1.0{PI/2.0*sinp.signum()}else{trig::asin(sinp)};
        let roll=trig::atan2(2.0*(q.w*q.x+q.y*q.z),1.0-2.0*(q.x*q.x+q.y*q.y));
        let yaw=trig::atan2(2.0*(q.w*q.z+q.x*q.y),1.0-2.0*(q.y*q.y+q.z*q.z));
        (roll,pitch,yaw)}
}

/// SLERP (Spherical Linear Interpolation).
pub fn slerp(a:Quaternion,b:Quaternion,t:f64)->Quaternion{
    let mut dot=a.dot(b);let mut b=b;
    if dot<0.0{b=b.scale(-1.0);dot=-dot;}
    if dot>0.9995{return a.add(b.add(a.scale(-1.0)).scale(t)).normalize();}
    let theta=trig::acos(dot);let sin_theta=trig::sin(theta);
    a.scale(trig::sin((1.0-t)*theta)/sin_theta).add(b.scale(trig::sin(t*theta)/sin_theta))
}

/// NLERP (Normalized Linear Interpolation, faster than SLERP).
pub fn nlerp(a:Quaternion,b:Quaternion,t:f64)->Quaternion{
    a.add(b.add(a.scale(-1.0)).scale(t)).normalize()
}

/// Squad (Spherical Quadrangle interpolation).
pub fn squad(q0:Quaternion,q1:Quaternion,s0:Quaternion,s1:Quaternion,t:f64)->Quaternion{
    slerp(slerp(q0,q1,t),slerp(s0,s1,t),2.0*t*(1.0-t))
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_identity(){let q=Quaternion::identity();let v=[1.0,2.0,3.0];let r=q.rotate_vector(&v);
        assert!(approx_eq(r[0],1.0,1e-10));assert!(approx_eq(r[1],2.0,1e-10));}
    #[test]fn test_rotation_90z(){let q=Quaternion::from_axis_angle(&[0.0,0.0,1.0],PI/2.0);
        let v=q.rotate_vector(&[1.0,0.0,0.0]);assert!(approx_eq(v[0],0.0,1e-10));assert!(approx_eq(v[1],1.0,1e-10));}
    #[test]fn test_norm(){let q=Quaternion::new(1.0,2.0,3.0,4.0).normalize();assert!(approx_eq(q.norm(),1.0,1e-10));}
    #[test]fn test_inverse(){let q=Quaternion::new(1.0,2.0,3.0,4.0);let r=q.mul(q.inverse());assert!(approx_eq(r.w/r.norm(),1.0,1e-6));}
    #[test]fn test_slerp_endpoints(){let a=Quaternion::identity();let b=Quaternion::from_axis_angle(&[0.0,0.0,1.0],PI);
        let s0=slerp(a,b,0.0);assert!(approx_eq(s0.w,a.w,1e-10));}
}
