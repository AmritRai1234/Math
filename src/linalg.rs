//! Linear algebra: vectors and matrices.

use crate::arithmetic::abs;
use crate::powers::sqrt;
use std::ops::{Add, Sub, Mul, Neg};

// ─── Vec2 ───
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vec2 { pub x: f64, pub y: f64 }

impl Vec2 {
    pub fn new(x: f64, y: f64) -> Self { Self { x, y } }
    pub fn zero() -> Self { Self { x: 0.0, y: 0.0 } }
    pub fn dot(self, o: Self) -> f64 { self.x*o.x + self.y*o.y }
    pub fn magnitude(self) -> f64 { sqrt(self.dot(self)) }
    pub fn normalize(self) -> Self { let m=self.magnitude(); Self{x:self.x/m,y:self.y/m} }
    pub fn scale(self, s: f64) -> Self { Self{x:self.x*s,y:self.y*s} }
    pub fn distance(self, o: Self) -> f64 { (self-o).magnitude() }
    pub fn angle_between(self, o: Self) -> f64 { crate::trig::acos(self.dot(o)/(self.magnitude()*o.magnitude())) }
    pub fn lerp(self, o: Self, t: f64) -> Self { Self{x:self.x+t*(o.x-self.x),y:self.y+t*(o.y-self.y)} }
    pub fn reflect(self, n: Self) -> Self { self - n.scale(2.0*self.dot(n)) }
    pub fn project(self, onto: Self) -> Self { onto.scale(self.dot(onto)/onto.dot(onto)) }
    pub fn reject(self, from: Self) -> Self { self - self.project(from) }
    pub fn rotate(self, angle: f64) -> Self {
        let (s,c) = crate::trig::sincos(angle);
        Self{x:self.x*c-self.y*s, y:self.x*s+self.y*c}
    }
}
impl Add for Vec2 { type Output=Self; fn add(self,o:Self)->Self{Self{x:self.x+o.x,y:self.y+o.y}} }
impl Sub for Vec2 { type Output=Self; fn sub(self,o:Self)->Self{Self{x:self.x-o.x,y:self.y-o.y}} }
impl Neg for Vec2 { type Output=Self; fn neg(self)->Self{Self{x:-self.x,y:-self.y}} }
impl Mul<f64> for Vec2 { type Output=Self; fn mul(self,s:f64)->Self{self.scale(s)} }

// ─── Vec3 ───
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vec3 { pub x: f64, pub y: f64, pub z: f64 }

impl Vec3 {
    pub fn new(x:f64,y:f64,z:f64)->Self{Self{x,y,z}}
    pub fn zero()->Self{Self{x:0.0,y:0.0,z:0.0}}
    pub fn dot(self,o:Self)->f64{self.x*o.x+self.y*o.y+self.z*o.z}
    pub fn cross(self,o:Self)->Self{Self{x:self.y*o.z-self.z*o.y,y:self.z*o.x-self.x*o.z,z:self.x*o.y-self.y*o.x}}
    pub fn magnitude(self)->f64{sqrt(self.dot(self))}
    pub fn normalize(self)->Self{let m=self.magnitude();Self{x:self.x/m,y:self.y/m,z:self.z/m}}
    pub fn scale(self,s:f64)->Self{Self{x:self.x*s,y:self.y*s,z:self.z*s}}
    pub fn distance(self,o:Self)->f64{(self-o).magnitude()}
    pub fn angle_between(self,o:Self)->f64{crate::trig::acos(self.dot(o)/(self.magnitude()*o.magnitude()))}
    pub fn lerp(self,o:Self,t:f64)->Self{Self{x:self.x+t*(o.x-self.x),y:self.y+t*(o.y-self.y),z:self.z+t*(o.z-self.z)}}
    pub fn reflect(self,n:Self)->Self{self-n.scale(2.0*self.dot(n))}
    pub fn refract(self,n:Self,eta:f64)->Self{let d=self.dot(n);let k=1.0-eta*eta*(1.0-d*d);if k<0.0{Self::zero()}else{self.scale(eta)-n.scale(eta*d+sqrt(k))}}
    pub fn project(self,onto:Self)->Self{onto.scale(self.dot(onto)/onto.dot(onto))}
    pub fn reject(self,from:Self)->Self{self-self.project(from)}
}
impl Add for Vec3{type Output=Self;fn add(self,o:Self)->Self{Self{x:self.x+o.x,y:self.y+o.y,z:self.z+o.z}}}
impl Sub for Vec3{type Output=Self;fn sub(self,o:Self)->Self{Self{x:self.x-o.x,y:self.y-o.y,z:self.z-o.z}}}
impl Neg for Vec3{type Output=Self;fn neg(self)->Self{Self{x:-self.x,y:-self.y,z:-self.z}}}
impl Mul<f64> for Vec3{type Output=Self;fn mul(self,s:f64)->Self{self.scale(s)}}

// ─── Vec4 ───
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vec4 { pub x:f64,pub y:f64,pub z:f64,pub w:f64 }
impl Vec4 {
    pub fn new(x:f64,y:f64,z:f64,w:f64)->Self{Self{x,y,z,w}}
    pub fn zero()->Self{Self{x:0.0,y:0.0,z:0.0,w:0.0}}
    pub fn dot(self,o:Self)->f64{self.x*o.x+self.y*o.y+self.z*o.z+self.w*o.w}
    pub fn magnitude(self)->f64{sqrt(self.dot(self))}
    pub fn normalize(self)->Self{let m=self.magnitude();Self{x:self.x/m,y:self.y/m,z:self.z/m,w:self.w/m}}
    pub fn scale(self,s:f64)->Self{Self{x:self.x*s,y:self.y*s,z:self.z*s,w:self.w*s}}
}
impl Add for Vec4{type Output=Self;fn add(self,o:Self)->Self{Self{x:self.x+o.x,y:self.y+o.y,z:self.z+o.z,w:self.w+o.w}}}
impl Sub for Vec4{type Output=Self;fn sub(self,o:Self)->Self{Self{x:self.x-o.x,y:self.y-o.y,z:self.z-o.z,w:self.w-o.w}}}
impl Neg for Vec4{type Output=Self;fn neg(self)->Self{Self{x:-self.x,y:-self.y,z:-self.z,w:-self.w}}}

// ─── Mat2 ───
#[derive(Clone,Copy,Debug)]
pub struct Mat2{pub m:[[f64;2];2]}
impl Mat2 {
    pub fn new(m:[[f64;2];2])->Self{Self{m}}
    pub fn identity()->Self{Self{m:[[1.0,0.0],[0.0,1.0]]}}
    pub fn zero()->Self{Self{m:[[0.0;2];2]}}
    pub fn determinant(self)->f64{self.m[0][0]*self.m[1][1]-self.m[0][1]*self.m[1][0]}
    pub fn transpose(self)->Self{Self{m:[[self.m[0][0],self.m[1][0]],[self.m[0][1],self.m[1][1]]]}}
    pub fn trace(self)->f64{self.m[0][0]+self.m[1][1]}
    pub fn inverse(self)->Option<Self>{let d=self.determinant();if abs(d)<1e-15{return None;}let inv_d=1.0/d;Some(Self{m:[[self.m[1][1]*inv_d,-self.m[0][1]*inv_d],[-self.m[1][0]*inv_d,self.m[0][0]*inv_d]]})}
    pub fn mul_vec(self,v:Vec2)->Vec2{Vec2{x:self.m[0][0]*v.x+self.m[0][1]*v.y,y:self.m[1][0]*v.x+self.m[1][1]*v.y}}
    pub fn scale(self,s:f64)->Self{Self{m:[[self.m[0][0]*s,self.m[0][1]*s],[self.m[1][0]*s,self.m[1][1]*s]]}}
}
impl Mul for Mat2{type Output=Self;fn mul(self,o:Self)->Self{let mut r=[[0.0;2];2];for i in 0..2{for j in 0..2{for k in 0..2{r[i][j]+=self.m[i][k]*o.m[k][j];}}}Self{m:r}}}

// ─── Mat3 ───
#[derive(Clone,Copy,Debug)]
pub struct Mat3{pub m:[[f64;3];3]}
impl Mat3 {
    pub fn new(m:[[f64;3];3])->Self{Self{m}}
    pub fn identity()->Self{let mut m=[[0.0;3];3];m[0][0]=1.0;m[1][1]=1.0;m[2][2]=1.0;Self{m}}
    pub fn zero()->Self{Self{m:[[0.0;3];3]}}
    pub fn transpose(self)->Self{let mut r=[[0.0;3];3];for i in 0..3{for j in 0..3{r[i][j]=self.m[j][i];}}Self{m:r}}
    pub fn trace(self)->f64{self.m[0][0]+self.m[1][1]+self.m[2][2]}
    pub fn determinant(self)->f64{
        self.m[0][0]*(self.m[1][1]*self.m[2][2]-self.m[1][2]*self.m[2][1])
       -self.m[0][1]*(self.m[1][0]*self.m[2][2]-self.m[1][2]*self.m[2][0])
       +self.m[0][2]*(self.m[1][0]*self.m[2][1]-self.m[1][1]*self.m[2][0])
    }
    pub fn inverse(self)->Option<Self>{let d=self.determinant();if abs(d)<1e-15{return None;}let id=1.0/d;let mut r=[[0.0;3];3];
        r[0][0]=(self.m[1][1]*self.m[2][2]-self.m[1][2]*self.m[2][1])*id;r[0][1]=(self.m[0][2]*self.m[2][1]-self.m[0][1]*self.m[2][2])*id;r[0][2]=(self.m[0][1]*self.m[1][2]-self.m[0][2]*self.m[1][1])*id;
        r[1][0]=(self.m[1][2]*self.m[2][0]-self.m[1][0]*self.m[2][2])*id;r[1][1]=(self.m[0][0]*self.m[2][2]-self.m[0][2]*self.m[2][0])*id;r[1][2]=(self.m[0][2]*self.m[1][0]-self.m[0][0]*self.m[1][2])*id;
        r[2][0]=(self.m[1][0]*self.m[2][1]-self.m[1][1]*self.m[2][0])*id;r[2][1]=(self.m[0][1]*self.m[2][0]-self.m[0][0]*self.m[2][1])*id;r[2][2]=(self.m[0][0]*self.m[1][1]-self.m[0][1]*self.m[1][0])*id;
        Some(Self{m:r})}
    pub fn mul_vec(self,v:Vec3)->Vec3{Vec3{x:self.m[0][0]*v.x+self.m[0][1]*v.y+self.m[0][2]*v.z,y:self.m[1][0]*v.x+self.m[1][1]*v.y+self.m[1][2]*v.z,z:self.m[2][0]*v.x+self.m[2][1]*v.y+self.m[2][2]*v.z}}
    pub fn scale(self,s:f64)->Self{let mut r=self.m;for i in 0..3{for j in 0..3{r[i][j]*=s;}}Self{m:r}}
}
impl Mul for Mat3{type Output=Self;fn mul(self,o:Self)->Self{let mut r=[[0.0;3];3];for i in 0..3{for j in 0..3{for k in 0..3{r[i][j]+=self.m[i][k]*o.m[k][j];}}}Self{m:r}}}

// ─── Mat4 ───
#[derive(Clone,Copy,Debug)]
pub struct Mat4{pub m:[[f64;4];4]}
impl Mat4 {
    pub fn new(m:[[f64;4];4])->Self{Self{m}}
    pub fn identity()->Self{let mut m=[[0.0;4];4];for i in 0..4{m[i][i]=1.0;}Self{m}}
    pub fn zero()->Self{Self{m:[[0.0;4];4]}}
    pub fn transpose(self)->Self{let mut r=[[0.0;4];4];for i in 0..4{for j in 0..4{r[i][j]=self.m[j][i];}}Self{m:r}}
    pub fn trace(self)->f64{self.m[0][0]+self.m[1][1]+self.m[2][2]+self.m[3][3]}
    pub fn mul_vec(self,v:Vec4)->Vec4{Vec4{x:self.m[0][0]*v.x+self.m[0][1]*v.y+self.m[0][2]*v.z+self.m[0][3]*v.w,y:self.m[1][0]*v.x+self.m[1][1]*v.y+self.m[1][2]*v.z+self.m[1][3]*v.w,z:self.m[2][0]*v.x+self.m[2][1]*v.y+self.m[2][2]*v.z+self.m[2][3]*v.w,w:self.m[3][0]*v.x+self.m[3][1]*v.y+self.m[3][2]*v.z+self.m[3][3]*v.w}}
    pub fn scale_mat(self,s:f64)->Self{let mut r=self.m;for i in 0..4{for j in 0..4{r[i][j]*=s;}}Self{m:r}}
    pub fn determinant(self)->f64{
        let m=&self.m;
        m[0][0]*(m[1][1]*(m[2][2]*m[3][3]-m[2][3]*m[3][2])-m[1][2]*(m[2][1]*m[3][3]-m[2][3]*m[3][1])+m[1][3]*(m[2][1]*m[3][2]-m[2][2]*m[3][1]))
       -m[0][1]*(m[1][0]*(m[2][2]*m[3][3]-m[2][3]*m[3][2])-m[1][2]*(m[2][0]*m[3][3]-m[2][3]*m[3][0])+m[1][3]*(m[2][0]*m[3][2]-m[2][2]*m[3][0]))
       +m[0][2]*(m[1][0]*(m[2][1]*m[3][3]-m[2][3]*m[3][1])-m[1][1]*(m[2][0]*m[3][3]-m[2][3]*m[3][0])+m[1][3]*(m[2][0]*m[3][1]-m[2][1]*m[3][0]))
       -m[0][3]*(m[1][0]*(m[2][1]*m[3][2]-m[2][2]*m[3][1])-m[1][1]*(m[2][0]*m[3][2]-m[2][2]*m[3][0])+m[1][2]*(m[2][0]*m[3][1]-m[2][1]*m[3][0]))
    }
}
impl Mul for Mat4{type Output=Self;fn mul(self,o:Self)->Self{let mut r=[[0.0;4];4];for i in 0..4{for j in 0..4{for k in 0..4{r[i][j]+=self.m[i][k]*o.m[k][j];}}}Self{m:r}}}

// ─── Transform helpers ───
pub fn rotate_2d(angle:f64)->Mat2{let(s,c)=crate::trig::sincos(angle);Mat2{m:[[c,-s],[s,c]]}}
pub fn scale_2d(sx:f64,sy:f64)->Mat2{Mat2{m:[[sx,0.0],[0.0,sy]]}}
pub fn rotate_x(a:f64)->Mat3{let(s,c)=crate::trig::sincos(a);Mat3{m:[[1.0,0.0,0.0],[0.0,c,-s],[0.0,s,c]]}}
pub fn rotate_y(a:f64)->Mat3{let(s,c)=crate::trig::sincos(a);Mat3{m:[[c,0.0,s],[0.0,1.0,0.0],[-s,0.0,c]]}}
pub fn rotate_z(a:f64)->Mat3{let(s,c)=crate::trig::sincos(a);Mat3{m:[[c,-s,0.0],[s,c,0.0],[0.0,0.0,1.0]]}}

pub fn translate_3d(tx:f64,ty:f64,tz:f64)->Mat4{let mut m=Mat4::identity();m.m[0][3]=tx;m.m[1][3]=ty;m.m[2][3]=tz;m}
pub fn scale_3d(sx:f64,sy:f64,sz:f64)->Mat4{let mut m=Mat4::identity();m.m[0][0]=sx;m.m[1][1]=sy;m.m[2][2]=sz;m}

pub fn look_at(eye:Vec3,target:Vec3,up:Vec3)->Mat4{
    let f=(target-eye).normalize(); let r=f.cross(up).normalize(); let u=r.cross(f);
    let mut m=Mat4::identity();
    m.m[0][0]=r.x;m.m[0][1]=r.y;m.m[0][2]=r.z;m.m[0][3]=-r.dot(eye);
    m.m[1][0]=u.x;m.m[1][1]=u.y;m.m[1][2]=u.z;m.m[1][3]=-u.dot(eye);
    m.m[2][0]=-f.x;m.m[2][1]=-f.y;m.m[2][2]=-f.z;m.m[2][3]=f.dot(eye);
    m
}

pub fn perspective(fov:f64,aspect:f64,near:f64,far:f64)->Mat4{
    let t=crate::trig::tan(fov/2.0);let mut m=Mat4::zero();
    m.m[0][0]=1.0/(aspect*t);m.m[1][1]=1.0/t;
    m.m[2][2]=-(far+near)/(far-near);m.m[2][3]=-2.0*far*near/(far-near);
    m.m[3][2]=-1.0;m
}

pub fn orthographic(l:f64,r:f64,b:f64,t:f64,n:f64,f:f64)->Mat4{
    let mut m=Mat4::identity();
    m.m[0][0]=2.0/(r-l);m.m[1][1]=2.0/(t-b);m.m[2][2]=-2.0/(f-n);
    m.m[0][3]=-(r+l)/(r-l);m.m[1][3]=-(t+b)/(t-b);m.m[2][3]=-(f+n)/(f-n);m
}

// ─── Decompositions (NxN dynamic) ───
pub fn lu_decompose(a:&[Vec<f64>])->(Vec<Vec<f64>>,Vec<Vec<f64>>){
    let n=a.len();let mut l=vec![vec![0.0;n];n];let mut u=vec![vec![0.0;n];n];
    for i in 0..n{l[i][i]=1.0;for j in i..n{let mut s=0.0;for k in 0..i{s+=l[i][k]*u[k][j];}u[i][j]=a[i][j]-s;}
    for j in(i+1)..n{let mut s=0.0;for k in 0..i{s+=l[j][k]*u[k][i];}l[j][i]=if abs(u[i][i])<1e-15{0.0}else{(a[j][i]-s)/u[i][i]};}}
    (l,u)
}

pub fn solve_linear(a:&[Vec<f64>],b:&[f64])->Vec<f64>{
    let n=a.len();let(l,u)=lu_decompose(a);
    let mut y=vec![0.0;n];for i in 0..n{let mut s=0.0;for j in 0..i{s+=l[i][j]*y[j];}y[i]=b[i]-s;}
    let mut x=vec![0.0;n];for i in(0..n).rev(){let mut s=0.0;for j in(i+1)..n{s+=u[i][j]*x[j];}x[i]=if abs(u[i][i])<1e-15{0.0}else{(y[i]-s)/u[i][i]};}
    x
}

pub fn outer_product(a:&[f64],b:&[f64])->Vec<Vec<f64>>{
    a.iter().map(|&ai|b.iter().map(|&bj|ai*bj).collect()).collect()
}

#[cfg(test)]
mod tests{
    use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_vec3_ops(){let a=Vec3::new(1.0,0.0,0.0);let b=Vec3::new(0.0,1.0,0.0);let c=a.cross(b);assert!(approx_eq(c.z,1.0,1e-14));assert!(approx_eq(a.dot(b),0.0,1e-14));}
    #[test]fn test_mat2_inv(){let m=Mat2::new([[1.0,2.0],[3.0,4.0]]);let inv=m.inverse().unwrap();let id=m*inv;assert!(approx_eq(id.m[0][0],1.0,1e-12));assert!(approx_eq(id.m[0][1],0.0,1e-12));}
    #[test]fn test_mat3_det(){let m=Mat3::new([[1.0,2.0,3.0],[0.0,1.0,4.0],[5.0,6.0,0.0]]);assert!(approx_eq(m.determinant(),1.0,1e-12));}
}
