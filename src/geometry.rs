//! 2D and 3D geometry: shapes, areas, volumes, distances, intersections.

use crate::arithmetic::abs;
use crate::constants;
use crate::powers;
use crate::trig;

// ─── 2D Shapes ───

pub fn circle_area(r:f64)->f64{constants::PI*r*r}
pub fn circle_circumference(r:f64)->f64{constants::TAU*r}
pub fn circle_arc_length(r:f64,angle:f64)->f64{r*angle}
pub fn circle_sector_area(r:f64,angle:f64)->f64{0.5*r*r*angle}

pub fn rect_area(w:f64,h:f64)->f64{w*h}
pub fn rect_perimeter(w:f64,h:f64)->f64{2.0*(w+h)}
pub fn rect_diagonal(w:f64,h:f64)->f64{powers::hypot(w,h)}

pub fn triangle_area_base_height(base:f64,height:f64)->f64{0.5*base*height}
pub fn triangle_area_sides(a:f64,b:f64,c:f64)->f64{
    let s=(a+b+c)/2.0;powers::sqrt(s*(s-a)*(s-b)*(s-c))// Heron's formula
}
pub fn triangle_area_vertices(x1:f64,y1:f64,x2:f64,y2:f64,x3:f64,y3:f64)->f64{
    abs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/2.0
}

pub fn regular_polygon_area(n:u32,side:f64)->f64{
    let nf=n as f64;nf*side*side/(4.0*trig::tan(constants::PI/nf))
}
pub fn regular_polygon_perimeter(n:u32,side:f64)->f64{n as f64*side}

pub fn ellipse_area(a:f64,b:f64)->f64{constants::PI*a*b}
pub fn ellipse_circumference_approx(a:f64,b:f64)->f64{
    // Ramanujan approximation
    let h=(a-b)*(a-b)/((a+b)*(a+b));
    constants::PI*(a+b)*(1.0+3.0*h/(10.0+powers::sqrt(4.0-3.0*h)))
}

pub fn trapezoid_area(a:f64,b:f64,h:f64)->f64{(a+b)*h/2.0}
pub fn parallelogram_area(base:f64,height:f64)->f64{base*height}
pub fn rhombus_area(d1:f64,d2:f64)->f64{d1*d2/2.0}

// ─── 3D Shapes ───

pub fn sphere_volume(r:f64)->f64{4.0/3.0*constants::PI*r*r*r}
pub fn sphere_surface_area(r:f64)->f64{4.0*constants::PI*r*r}

pub fn cylinder_volume(r:f64,h:f64)->f64{constants::PI*r*r*h}
pub fn cylinder_surface_area(r:f64,h:f64)->f64{constants::TAU*r*(r+h)}
pub fn cylinder_lateral_area(r:f64,h:f64)->f64{constants::TAU*r*h}

pub fn cone_volume(r:f64,h:f64)->f64{constants::PI*r*r*h/3.0}
pub fn cone_surface_area(r:f64,h:f64)->f64{let s=powers::sqrt(r*r+h*h);constants::PI*r*(r+s)}
pub fn cone_slant_height(r:f64,h:f64)->f64{powers::sqrt(r*r+h*h)}

pub fn cube_volume(s:f64)->f64{s*s*s}
pub fn cube_surface_area(s:f64)->f64{6.0*s*s}
pub fn cube_diagonal(s:f64)->f64{s*powers::sqrt(3.0)}

pub fn box_volume(l:f64,w:f64,h:f64)->f64{l*w*h}
pub fn box_surface_area(l:f64,w:f64,h:f64)->f64{2.0*(l*w+w*h+h*l)}
pub fn box_diagonal(l:f64,w:f64,h:f64)->f64{powers::sqrt(l*l+w*w+h*h)}

pub fn torus_volume(major_r:f64,minor_r:f64)->f64{constants::TAU*constants::PI*major_r*minor_r*minor_r}
pub fn torus_surface_area(major_r:f64,minor_r:f64)->f64{constants::TAU*constants::TAU*major_r*minor_r}

pub fn pyramid_volume(base_area:f64,h:f64)->f64{base_area*h/3.0}
pub fn frustum_volume(r1:f64,r2:f64,h:f64)->f64{constants::PI*h*(r1*r1+r2*r2+r1*r2)/3.0}

pub fn ellipsoid_volume(a:f64,b:f64,c:f64)->f64{4.0/3.0*constants::PI*a*b*c}

// ─── Distance & Intersection ───

pub fn point_distance_2d(x1:f64,y1:f64,x2:f64,y2:f64)->f64{powers::hypot(x2-x1,y2-y1)}
pub fn point_distance_3d(x1:f64,y1:f64,z1:f64,x2:f64,y2:f64,z2:f64)->f64{powers::sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))}
pub fn manhattan_distance(x1:f64,y1:f64,x2:f64,y2:f64)->f64{abs(x2-x1)+abs(y2-y1)}
pub fn chebyshev_distance(x1:f64,y1:f64,x2:f64,y2:f64)->f64{abs(x2-x1).max(abs(y2-y1))}
pub fn minkowski_distance(x1:f64,y1:f64,x2:f64,y2:f64,p:f64)->f64{(abs(x2-x1).powf(p)+abs(y2-y1).powf(p)).powf(1.0/p)}

pub fn point_to_line_distance(px:f64,py:f64,x1:f64,y1:f64,x2:f64,y2:f64)->f64{
    abs((y2-y1)*px-(x2-x1)*py+x2*y1-y2*x1)/point_distance_2d(x1,y1,x2,y2)
}

/// Check if point is inside triangle (using barycentric coordinates).
pub fn point_in_triangle(px:f64,py:f64,x1:f64,y1:f64,x2:f64,y2:f64,x3:f64,y3:f64)->bool{
    let d=(y2-y3)*(x1-x3)+(x3-x2)*(y1-y3);if abs(d)<1e-15{return false;}
    let a=((y2-y3)*(px-x3)+(x3-x2)*(py-y3))/d;
    let b=((y3-y1)*(px-x3)+(x1-x3)*(py-y3))/d;
    let c=1.0-a-b;
    a>=0.0&&b>=0.0&&c>=0.0
}

/// Check if point is inside circle.
pub fn point_in_circle(px:f64,py:f64,cx:f64,cy:f64,r:f64)->bool{
    (px-cx)*(px-cx)+(py-cy)*(py-cy)<=r*r
}

/// Check if two circles intersect.
pub fn circles_intersect(x1:f64,y1:f64,r1:f64,x2:f64,y2:f64,r2:f64)->bool{
    let d=point_distance_2d(x1,y1,x2,y2);
    d<=r1+r2&&d>=abs(r1-r2)
}

/// Line-line intersection. Returns Some((x,y)) or None if parallel.
pub fn line_intersection(x1:f64,y1:f64,x2:f64,y2:f64,x3:f64,y3:f64,x4:f64,y4:f64)->Option<(f64,f64)>{
    let d=(x1-x2)*(y3-y4)-(y1-y2)*(x3-x4);
    if abs(d)<1e-15{return None;}
    let t=((x1-x3)*(y3-y4)-(y1-y3)*(x3-x4))/d;
    Some((x1+t*(x2-x1),y1+t*(y2-y1)))
}

/// Convex hull (Andrew's monotone chain). Returns vertices in CCW order.
pub fn convex_hull(points:&[(f64,f64)])->Vec<(f64,f64)>{
    let mut pts=points.to_vec();
    pts.sort_by(|a,b|a.0.partial_cmp(&b.0).unwrap().then(a.1.partial_cmp(&b.1).unwrap()));
    let n=pts.len();if n<=1{return pts;}
    let mut hull=Vec::with_capacity(2*n);
    for &p in &pts{while hull.len()>=2&&cross_2d(hull[hull.len()-2],hull[hull.len()-1],p)<=0.0{hull.pop();}hull.push(p);}
    let lower=hull.len()+1;
    for &p in pts.iter().rev(){while hull.len()>=lower&&cross_2d(hull[hull.len()-2],hull[hull.len()-1],p)<=0.0{hull.pop();}hull.push(p);}
    hull.pop();hull
}

fn cross_2d(o:(f64,f64),a:(f64,f64),b:(f64,f64))->f64{(a.0-o.0)*(b.1-o.1)-(a.1-o.1)*(b.0-o.0)}

/// Polygon area (shoelace formula).
pub fn polygon_area(vertices:&[(f64,f64)])->f64{
    let n=vertices.len();if n<3{return 0.0;}
    let mut sum=0.0;for i in 0..n{let j=(i+1)%n;sum+=vertices[i].0*vertices[j].1-vertices[j].0*vertices[i].1;}abs(sum)/2.0
}

/// Centroid of polygon.
pub fn polygon_centroid(vertices:&[(f64,f64)])->(f64,f64){
    let n=vertices.len()as f64;let sx:f64=vertices.iter().map(|v|v.0).sum();let sy:f64=vertices.iter().map(|v|v.1).sum();(sx/n,sy/n)
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_circle(){assert!(approx_eq(circle_area(1.0),constants::PI,1e-12));}
    #[test]fn test_sphere(){assert!(approx_eq(sphere_volume(1.0),4.0/3.0*constants::PI,1e-12));}
    #[test]fn test_heron(){assert!(approx_eq(triangle_area_sides(3.0,4.0,5.0),6.0,1e-10));}
    #[test]fn test_point_in_triangle(){assert!(point_in_triangle(0.5,0.5,0.0,0.0,2.0,0.0,0.0,2.0));assert!(!point_in_triangle(3.0,3.0,0.0,0.0,2.0,0.0,0.0,2.0));}
    #[test]fn test_polygon_area(){let sq=[(0.0,0.0),(1.0,0.0),(1.0,1.0),(0.0,1.0)];assert!(approx_eq(polygon_area(&sq),1.0,1e-12));}
    #[test]fn test_convex_hull(){let pts=[(0.0,0.0),(1.0,0.0),(0.5,0.5),(0.0,1.0),(1.0,1.0)];let hull=convex_hull(&pts);assert_eq!(hull.len(),4);}
}
