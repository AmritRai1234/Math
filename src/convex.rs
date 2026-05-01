//! Convex geometry: convex hull, Minkowski sum, point-in-polygon.

/// Convex hull (2D, Graham scan). Returns indices of hull vertices in CCW order.
pub fn convex_hull(points:&[(f64,f64)])->Vec<usize>{
    let n=points.len();if n<3{return(0..n).collect();}
    // Find lowest point
    let mut pivot=0;for i in 1..n{if points[i].1<points[pivot].1||(points[i].1==points[pivot].1&&points[i].0<points[pivot].0){pivot=i;}}
    let mut indices:Vec<usize>=(0..n).filter(|&i|i!=pivot).collect();
    let p0=points[pivot];
    indices.sort_by(|&a,&b|{let ang_a=crate::trig::atan2(points[a].1-p0.1,points[a].0-p0.0);
        let ang_b=crate::trig::atan2(points[b].1-p0.1,points[b].0-p0.0);ang_a.partial_cmp(&ang_b).unwrap()});
    let mut hull=vec![pivot,indices[0]];
    for i in 1..indices.len(){
        while hull.len()>1{let(a,b,c)=(hull[hull.len()-2],hull[hull.len()-1],indices[i]);
            let cross=(points[b].0-points[a].0)*(points[c].1-points[a].1)-(points[b].1-points[a].1)*(points[c].0-points[a].0);
            if cross<=0.0{hull.pop();}else{break;}}
        hull.push(indices[i]);}hull
}

/// Polygon area (signed, CCW positive).
pub fn polygon_area(vertices:&[(f64,f64)])->f64{
    let n=vertices.len();let mut area=0.0;
    for i in 0..n{let j=(i+1)%n;area+=vertices[i].0*vertices[j].1-vertices[j].0*vertices[i].1;}area/2.0
}

/// Point in convex polygon test (using cross products).
pub fn point_in_polygon(point:(f64,f64),polygon:&[(f64,f64)])->bool{
    let n=polygon.len();let mut pos=0;let mut neg=0;
    for i in 0..n{let j=(i+1)%n;
        let cross=(polygon[j].0-polygon[i].0)*(point.1-polygon[i].1)-(polygon[j].1-polygon[i].1)*(point.0-polygon[i].0);
        if cross>0.0{pos+=1;}else if cross<0.0{neg+=1;}
        if pos>0&&neg>0{return false;}}true
}

/// Minkowski sum of two convex polygons.
pub fn minkowski_sum(a:&[(f64,f64)],b:&[(f64,f64)])->Vec<(f64,f64)>{
    let mut result=Vec::new();
    for &(ax,ay) in a{for &(bx,by) in b{result.push((ax+bx,ay+by));}}
    let indices=convex_hull(&result);indices.iter().map(|&i|result[i]).collect()
}

/// Centroid of polygon.
pub fn polygon_centroid(vertices:&[(f64,f64)])->(f64,f64){
    let n=vertices.len();let area=polygon_area(vertices);let mut cx=0.0;let mut cy=0.0;
    for i in 0..n{let j=(i+1)%n;let cross=vertices[i].0*vertices[j].1-vertices[j].0*vertices[i].1;
        cx+=(vertices[i].0+vertices[j].0)*cross;cy+=(vertices[i].1+vertices[j].1)*cross;}
    (cx/(6.0*area),cy/(6.0*area))
}

/// Perimeter of polygon.
pub fn polygon_perimeter(vertices:&[(f64,f64)])->f64{
    let n=vertices.len();(0..n).map(|i|{let j=(i+1)%n;
        let(dx,dy)=(vertices[j].0-vertices[i].0,vertices[j].1-vertices[i].1);crate::powers::sqrt(dx*dx+dy*dy)}).sum()
}

/// Distance from point to line segment.
pub fn point_to_segment_distance(p:(f64,f64),a:(f64,f64),b:(f64,f64))->f64{
    let(dx,dy)=(b.0-a.0,b.1-a.1);let len2=dx*dx+dy*dy;
    if len2<1e-30{return crate::powers::sqrt((p.0-a.0)*(p.0-a.0)+(p.1-a.1)*(p.1-a.1));}
    let t=((p.0-a.0)*dx+(p.1-a.1)*dy)/len2;let t=t.max(0.0).min(1.0);
    let proj=(a.0+t*dx,a.1+t*dy);crate::powers::sqrt((p.0-proj.0)*(p.0-proj.0)+(p.1-proj.1)*(p.1-proj.1))
}

/// Line segment intersection test.
pub fn segments_intersect(a1:(f64,f64),a2:(f64,f64),b1:(f64,f64),b2:(f64,f64))->bool{
    let d1=cross_2d(b1,b2,a1);let d2=cross_2d(b1,b2,a2);let d3=cross_2d(a1,a2,b1);let d4=cross_2d(a1,a2,b2);
    if((d1>0.0&&d2<0.0)||(d1<0.0&&d2>0.0))&&((d3>0.0&&d4<0.0)||(d3<0.0&&d4>0.0)){return true;}false
}
fn cross_2d(a:(f64,f64),b:(f64,f64),c:(f64,f64))->f64{(b.0-a.0)*(c.1-a.1)-(b.1-a.1)*(c.0-a.0)}

/// Closest pair of points (brute force).
pub fn closest_pair(points:&[(f64,f64)])->((f64,f64),(f64,f64),f64){
    let mut min_d=f64::INFINITY;let mut best=((0.0,0.0),(0.0,0.0));
    for i in 0..points.len(){for j in i+1..points.len(){
        let d=crate::powers::sqrt((points[i].0-points[j].0)*(points[i].0-points[j].0)+(points[i].1-points[j].1)*(points[i].1-points[j].1));
        if d<min_d{min_d=d;best=(points[i],points[j]);}}}(best.0,best.1,min_d)
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_hull(){let pts=vec![(0.0,0.0),(1.0,0.0),(0.5,0.5),(0.0,1.0),(1.0,1.0)];
        let hull=convex_hull(&pts);assert_eq!(hull.len(),4);}
    #[test]fn test_area(){let sq=vec![(0.0,0.0),(1.0,0.0),(1.0,1.0),(0.0,1.0)];
        assert!(approx_eq(polygon_area(&sq),1.0,1e-10));}
    #[test]fn test_pip(){let sq=vec![(0.0,0.0),(1.0,0.0),(1.0,1.0),(0.0,1.0)];
        assert!(point_in_polygon((0.5,0.5),&sq));}
}
