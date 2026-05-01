//! Calendar and date mathematics.

/// Julian Day Number from (year, month, day).
pub fn julian_day(year:i32,month:u32,day:u32)->i64{
    let a=(14-month as i32)/12;let y=year+4800-a;let m=month as i32+12*a-3;
    day as i64+((153*m+2)/5)as i64+365*y as i64+(y/4)as i64-(y/100)as i64+(y/400)as i64-32045
}

/// Julian Day to (year, month, day).
pub fn from_julian_day(jd:i64)->(i32,u32,u32){
    let a=jd+32044;let b=(4*a+3)/146097;let c=a-(146097*b)/4;
    let d=(4*c+3)/1461;let e=c-(1461*d)/4;let m=(5*e+2)/153;
    let day=(e-(153*m+2)/5+1)as u32;let month=(m+3-12*(m/10))as u32;
    let year=(100*b+d-4800+m/10)as i32;(year,month,day)
}

/// Day of week (0=Monday, 6=Sunday) using Zeller's congruence.
pub fn day_of_week(year:i32,month:u32,day:u32)->u32{
    let jd=julian_day(year,month,day);(jd as u32 %7) // 0=Mon for this JD formula
}

/// Day of week name.
pub fn day_name(year:i32,month:u32,day:u32)->&'static str{
    ["Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday"][day_of_week(year,month,day)as usize]
}

/// Is leap year (Gregorian).
pub fn is_leap_year(year:i32)->bool{year%4==0&&(year%100!=0||year%400==0)}

/// Days in month.
pub fn days_in_month(year:i32,month:u32)->u32{
    match month{1|3|5|7|8|10|12=>31,4|6|9|11=>30,2=>if is_leap_year(year){29}else{28},_=>0}
}

/// Days in year.
pub fn days_in_year(year:i32)->u32{if is_leap_year(year){366}else{365}}

/// Day of year (1-based).
pub fn day_of_year(year:i32,month:u32,day:u32)->u32{
    (1..month).map(|m|days_in_month(year,m)).sum::<u32>()+day
}

/// ISO week number.
pub fn iso_week(year:i32,month:u32,day:u32)->u32{
    let doy=day_of_year(year,month,day);let dow=day_of_week(year,month,day);
    let week=((doy+6)as i32-dow as i32)/7+1;if week<1{52}else{week as u32}
}

/// Days between two dates.
pub fn days_between(y1:i32,m1:u32,d1:u32,y2:i32,m2:u32,d2:u32)->i64{
    julian_day(y2,m2,d2)-julian_day(y1,m1,d1)
}

/// Easter date (Computus — Anonymous Gregorian algorithm).
pub fn easter(year:i32)->(u32,u32){
    let a=year%19;let b=year/100;let c=year%100;let d=b/4;let e=b%4;
    let f=(b+8)/25;let g=(b-f+1)/3;let h=(19*a+b-d-g+15)%30;
    let i=c/4;let k=c%4;let l=(32+2*e+2*i-h-k)%7;let m=(a+11*h+22*l)/451;
    let month=(h+l-7*m+114)/31;let day=(h+l-7*m+114)%31+1;
    (month as u32,day as u32)
}

/// Unix timestamp to (year, month, day, hour, minute, second).
pub fn from_unix(timestamp:i64)->(i32,u32,u32,u32,u32,u32){
    let days=timestamp/86400;let time_of_day=(timestamp%86400+86400)%86400;
    let jd=days+2440588; // Unix epoch = JD 2440587.5
    let(y,m,d)=from_julian_day(jd);
    let h=(time_of_day/3600)as u32;let min=((time_of_day%3600)/60)as u32;let s=(time_of_day%60)as u32;
    (y,m,d,h,min,s)
}

/// (year, month, day) to Unix timestamp (midnight UTC).
pub fn to_unix(year:i32,month:u32,day:u32)->i64{
    (julian_day(year,month,day)-2440588)*86400
}

/// Age in years between two dates.
pub fn age_in_years(by:i32,bm:u32,bd:u32,cy:i32,cm:u32,cd:u32)->u32{
    let mut age=(cy-by)as u32;if cm<bm||(cm==bm&&cd<bd){age-=1;}age
}

/// Modified Julian Day.
pub fn modified_julian_day(year:i32,month:u32,day:u32)->f64{julian_day(year,month,day)as f64-2400000.5}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_julian(){assert_eq!(julian_day(2000,1,1),2451545);}
    #[test]fn test_roundtrip(){let(y,m,d)=from_julian_day(julian_day(2024,7,15));assert_eq!((y,m,d),(2024,7,15));}
    #[test]fn test_dow(){assert_eq!(day_name(2024,1,1),"Monday");}
    #[test]fn test_leap(){assert!(is_leap_year(2024));assert!(!is_leap_year(1900));assert!(is_leap_year(2000));}
    #[test]fn test_easter(){assert_eq!(easter(2024),(3,31));} // March 31
    #[test]fn test_days_between(){assert_eq!(days_between(2024,1,1,2024,12,31),365);}
}
