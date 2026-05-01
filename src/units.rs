//! Unit conversions and dimensional analysis.

// ─── Length ───
pub fn meters_to_feet(m:f64)->f64{m*3.28084}
pub fn feet_to_meters(ft:f64)->f64{ft/3.28084}
pub fn km_to_miles(km:f64)->f64{km*0.621371}
pub fn miles_to_km(mi:f64)->f64{mi/0.621371}
pub fn inches_to_cm(inch:f64)->f64{inch*2.54}
pub fn cm_to_inches(cm:f64)->f64{cm/2.54}
pub fn nautical_miles_to_km(nm:f64)->f64{nm*1.852}
pub fn light_years_to_km(ly:f64)->f64{ly*9.461e12}

// ─── Mass ───
pub fn kg_to_lbs(kg:f64)->f64{kg*2.20462}
pub fn lbs_to_kg(lbs:f64)->f64{lbs/2.20462}
pub fn kg_to_oz(kg:f64)->f64{kg*35.274}
pub fn grams_to_grains(g:f64)->f64{g*15.4324}

// ─── Temperature ───
pub fn celsius_to_fahrenheit(c:f64)->f64{c*9.0/5.0+32.0}
pub fn fahrenheit_to_celsius(f:f64)->f64{(f-32.0)*5.0/9.0}
pub fn celsius_to_kelvin(c:f64)->f64{c+273.15}
pub fn kelvin_to_celsius(k:f64)->f64{k-273.15}
pub fn fahrenheit_to_kelvin(f:f64)->f64{celsius_to_kelvin(fahrenheit_to_celsius(f))}

// ─── Speed ───
pub fn kmh_to_ms(kmh:f64)->f64{kmh/3.6}
pub fn ms_to_kmh(ms:f64)->f64{ms*3.6}
pub fn mph_to_kmh(mph:f64)->f64{mph*1.60934}
pub fn knots_to_kmh(kn:f64)->f64{kn*1.852}

// ─── Area ───
pub fn sqm_to_sqft(sqm:f64)->f64{sqm*10.7639}
pub fn acres_to_sqm(acres:f64)->f64{acres*4046.86}
pub fn hectares_to_acres(ha:f64)->f64{ha*2.47105}

// ─── Volume ───
pub fn liters_to_gallons(l:f64)->f64{l*0.264172}
pub fn gallons_to_liters(gal:f64)->f64{gal/0.264172}
pub fn ml_to_fl_oz(ml:f64)->f64{ml*0.033814}
pub fn cubic_m_to_liters(m3:f64)->f64{m3*1000.0}

// ─── Pressure ───
pub fn atm_to_pascal(atm:f64)->f64{atm*101325.0}
pub fn pascal_to_atm(pa:f64)->f64{pa/101325.0}
pub fn psi_to_pascal(psi:f64)->f64{psi*6894.76}
pub fn bar_to_pascal(bar:f64)->f64{bar*100000.0}

// ─── Energy ───
pub fn joules_to_calories(j:f64)->f64{j/4.184}
pub fn calories_to_joules(cal:f64)->f64{cal*4.184}
pub fn kwh_to_joules(kwh:f64)->f64{kwh*3.6e6}
pub fn ev_to_joules(ev:f64)->f64{ev*1.602176634e-19}

// ─── Angle ───
pub fn degrees_to_radians(deg:f64)->f64{deg*crate::constants::PI/180.0}
pub fn radians_to_degrees(rad:f64)->f64{rad*180.0/crate::constants::PI}
pub fn degrees_to_gradians(deg:f64)->f64{deg*10.0/9.0}
pub fn gradians_to_degrees(grad:f64)->f64{grad*9.0/10.0}
pub fn turns_to_radians(turns:f64)->f64{turns*crate::constants::TAU}

// ─── Time ───
pub fn hours_to_seconds(h:f64)->f64{h*3600.0}
pub fn days_to_seconds(d:f64)->f64{d*86400.0}
pub fn years_to_seconds(y:f64)->f64{y*365.25*86400.0}

// ─── Data ───
pub fn bytes_to_bits(b:f64)->f64{b*8.0}
pub fn kb_to_bytes(kb:f64)->f64{kb*1024.0}
pub fn mb_to_bytes(mb:f64)->f64{mb*1048576.0}
pub fn gb_to_bytes(gb:f64)->f64{gb*1073741824.0}

// ─── Physical Constants ───
pub const SPEED_OF_LIGHT:f64=299792458.0;      // m/s
pub const PLANCK:f64=6.62607015e-34;            // J·s
pub const BOLTZMANN:f64=1.380649e-23;           // J/K
pub const AVOGADRO:f64=6.02214076e23;           // 1/mol
pub const GRAVITATIONAL:f64=6.67430e-11;        // m³/(kg·s²)
pub const ELECTRON_MASS:f64=9.1093837015e-31;   // kg
pub const PROTON_MASS:f64=1.67262192369e-27;    // kg
pub const ELEMENTARY_CHARGE:f64=1.602176634e-19; // C
pub const GAS_CONSTANT:f64=8.314462618;         // J/(mol·K)

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_temp(){assert!(approx_eq(celsius_to_fahrenheit(100.0),212.0,1e-10));assert!(approx_eq(fahrenheit_to_celsius(32.0),0.0,1e-10));}
    #[test]fn test_length(){assert!(approx_eq(km_to_miles(1.0),0.621371,1e-4));}
    #[test]fn test_roundtrip(){let m=100.0;assert!(approx_eq(feet_to_meters(meters_to_feet(m)),m,1e-6));}
}
