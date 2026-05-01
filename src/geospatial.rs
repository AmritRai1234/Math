//! Map projections and geospatial math.

use crate::trig;use crate::exponential;use crate::powers;use crate::constants::PI;

/// Mercator projection: (lat, lon) → (x, y). Input in radians.
pub fn mercator(lat:f64,lon:f64)->(f64,f64){(lon,exponential::ln(trig::tan(PI/4.0+lat/2.0)))}
/// Inverse Mercator.
pub fn mercator_inv(x:f64,y:f64)->(f64,f64){(2.0*trig::atan(exponential::exp(y))-PI/2.0,x)}

/// Web Mercator (EPSG:3857). Input degrees, output meters.
pub fn web_mercator(lat_deg:f64,lon_deg:f64)->(f64,f64){
    let r=6378137.0;let lon=lon_deg*PI/180.0;let lat=lat_deg*PI/180.0;
    (r*lon,r*exponential::ln(trig::tan(PI/4.0+lat/2.0)))
}

/// Stereographic projection.
pub fn stereographic(lat:f64,lon:f64)->(f64,f64){
    let k=2.0/(1.0+trig::cos(lat)*trig::cos(lon));(k*trig::cos(lat)*trig::sin(lon),k*trig::sin(lat))
}

/// Lambert conformal conic (simplified, single parallel).
pub fn lambert_conformal(lat:f64,lon:f64,lat0:f64,lon0:f64)->(f64,f64){
    let n=trig::sin(lat0);let rho=1.0/trig::tan(PI/4.0+lat/2.0);
    let rho0=1.0/trig::tan(PI/4.0+lat0/2.0);
    let theta=n*(lon-lon0);
    (rho*trig::sin(theta),rho0-rho*trig::cos(theta))
}

/// Equirectangular projection (plate carrée).
pub fn equirectangular(lat:f64,lon:f64)->(f64,f64){(lon,lat)}

/// Azimuthal equidistant projection.
pub fn azimuthal_equidistant(lat:f64,lon:f64,lat0:f64,lon0:f64)->(f64,f64){
    let c=trig::acos(trig::sin(lat0)*trig::sin(lat)+trig::cos(lat0)*trig::cos(lat)*trig::cos(lon-lon0));
    if c.abs()<1e-15{return(0.0,0.0);}
    let k=c/trig::sin(c);
    (k*trig::cos(lat)*trig::sin(lon-lon0),k*(trig::cos(lat0)*trig::sin(lat)-trig::sin(lat0)*trig::cos(lat)*trig::cos(lon-lon0)))
}

/// Haversine distance (meters) between two points (lat/lon in degrees).
pub fn haversine_distance(lat1:f64,lon1:f64,lat2:f64,lon2:f64)->f64{
    let r=6371000.0;let(la1,lo1,la2,lo2)=(lat1*PI/180.0,lon1*PI/180.0,lat2*PI/180.0,lon2*PI/180.0);
    let dlat=la2-la1;let dlon=lo2-lo1;
    let a=trig::sin(dlat/2.0)*trig::sin(dlat/2.0)+trig::cos(la1)*trig::cos(la2)*trig::sin(dlon/2.0)*trig::sin(dlon/2.0);
    r*2.0*trig::atan2(powers::sqrt(a),powers::sqrt(1.0-a))
}

/// Bearing from point 1 to point 2 (degrees).
pub fn bearing(lat1:f64,lon1:f64,lat2:f64,lon2:f64)->f64{
    let(la1,lo1,la2,lo2)=(lat1*PI/180.0,lon1*PI/180.0,lat2*PI/180.0,lon2*PI/180.0);
    let dlon=lo2-lo1;let y=trig::sin(dlon)*trig::cos(la2);
    let x=trig::cos(la1)*trig::sin(la2)-trig::sin(la1)*trig::cos(la2)*trig::cos(dlon);
    (trig::atan2(y,x)*180.0/PI+360.0)%360.0
}

/// Destination point given start, bearing (degrees), and distance (meters).
pub fn destination_point(lat:f64,lon:f64,bearing_deg:f64,distance:f64)->(f64,f64){
    let r=6371000.0;let la=lat*PI/180.0;let lo=lon*PI/180.0;let brng=bearing_deg*PI/180.0;let d=distance/r;
    let lat2=trig::asin(trig::sin(la)*trig::cos(d)+trig::cos(la)*trig::sin(d)*trig::cos(brng));
    let lon2=lo+trig::atan2(trig::sin(brng)*trig::sin(d)*trig::cos(la),trig::cos(d)-trig::sin(la)*trig::sin(lat2));
    (lat2*180.0/PI,lon2*180.0/PI)
}

/// Midpoint between two coordinates (degrees).
pub fn midpoint(lat1:f64,lon1:f64,lat2:f64,lon2:f64)->(f64,f64){
    let(la1,lo1,la2,lo2)=(lat1*PI/180.0,lon1*PI/180.0,lat2*PI/180.0,lon2*PI/180.0);
    let bx=trig::cos(la2)*trig::cos(lo2-lo1);let by=trig::cos(la2)*trig::sin(lo2-lo1);
    let lat3=trig::atan2(trig::sin(la1)+trig::sin(la2),powers::sqrt((trig::cos(la1)+bx)*(trig::cos(la1)+bx)+by*by));
    let lon3=lo1+trig::atan2(by,trig::cos(la1)+bx);
    (lat3*180.0/PI,lon3*180.0/PI)
}

/// UTM zone from longitude.
pub fn utm_zone(lon_deg:f64)->u32{((lon_deg+180.0)/6.0).floor()as u32+1}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_haversine(){let d=haversine_distance(51.5074,-0.1278,48.8566,2.3522);assert!(d>340000.0&&d<345000.0);}
    #[test]fn test_mercator_roundtrip(){let(x,y)=mercator(0.5,1.0);let(lat,lon)=mercator_inv(x,y);
        assert!(approx_eq(lat,0.5,1e-6));assert!(approx_eq(lon,1.0,1e-6));}
    #[test]fn test_utm(){assert_eq!(utm_zone(-73.9857),18);} // New York
}
