//! Time series analysis.

use crate::arithmetic::abs;
use crate::stats;

/// Simple moving average.
pub fn sma(data:&[f64],window:usize)->Vec<f64>{
    if data.len()<window{return vec![];}
    let mut result=Vec::with_capacity(data.len()-window+1);let mut sum:f64=data[..window].iter().sum();
    result.push(sum/window as f64);
    for i in window..data.len(){sum+=data[i]-data[i-window];result.push(sum/window as f64);}result
}

/// Exponential moving average.
pub fn ema(data:&[f64],alpha:f64)->Vec<f64>{
    if data.is_empty(){return vec![];}
    let mut result=vec![data[0]];
    for i in 1..data.len(){result.push(alpha*data[i]+(1.0-alpha)*result[i-1]);}result
}

/// Double exponential smoothing (Holt's method).
pub fn double_exp_smooth(data:&[f64],alpha:f64,beta:f64)->Vec<f64>{
    if data.len()<2{return data.to_vec();}
    let mut level=data[0];let mut trend=data[1]-data[0];let mut result=vec![level];
    for i in 1..data.len(){let prev_level=level;level=alpha*data[i]+(1.0-alpha)*(prev_level+trend);
        trend=beta*(level-prev_level)+(1.0-beta)*trend;result.push(level+trend);}result
}

/// Autocorrelation at lag k.
pub fn autocorrelation(data:&[f64],lag:usize)->f64{
    let n=data.len();if lag>=n{return 0.0;}let mean=stats::mean(data);
    let var:f64=data.iter().map(|&x|(x-mean)*(x-mean)).sum();if var==0.0{return 0.0;}
    let cov:f64=(0..n-lag).map(|i|(data[i]-mean)*(data[i+lag]-mean)).sum();cov/var
}

/// Autocorrelation function for lags 0..max_lag.
pub fn acf(data:&[f64],max_lag:usize)->Vec<f64>{(0..=max_lag).map(|k|autocorrelation(data,k)).collect()}

/// Partial autocorrelation (Durbin-Levinson).
pub fn pacf(data:&[f64],max_lag:usize)->Vec<f64>{
    let r:Vec<f64>=acf(data,max_lag);let mut phi=vec![vec![0.0;max_lag+1];max_lag+1];
    let mut result=vec![1.0];
    if max_lag==0{return result;}
    phi[1][1]=r[1];result.push(r[1]);
    for k in 2..=max_lag{let mut num=r[k];let mut den=1.0;
        for j in 1..k{num-=phi[k-1][j]*r[k-j];den-=phi[k-1][j]*r[j];}
        if abs(den)<1e-30{result.push(0.0);continue;}
        phi[k][k]=num/den;result.push(phi[k][k]);
        for j in 1..k{phi[k][j]=phi[k-1][j]-phi[k][k]*phi[k-1][k-j];}}
    result
}

/// Differencing (d times).
pub fn difference(data:&[f64],d:usize)->Vec<f64>{
    let mut result=data.to_vec();
    for _ in 0..d{let prev=result.clone();result=prev.windows(2).map(|w|w[1]-w[0]).collect();}result
}

/// Augmented Dickey-Fuller test statistic (simplified).
pub fn adf_statistic(data:&[f64])->f64{
    let diff=difference(data,1);let n=diff.len();if n<2{return 0.0;}
    let y:Vec<f64>=diff.iter().cloned().collect();let x:Vec<f64>=data[..n].to_vec();
    let(a,b)=crate::regression::linear_regression(&x,&y);
    let residuals:Vec<f64>=x.iter().zip(y.iter()).map(|(&xi,&yi)|yi-(a+b*xi)).collect();
    let se=crate::powers::sqrt(residuals.iter().map(|&r|r*r).sum::<f64>()/(n as f64-2.0));
    let x_var:f64=x.iter().map(|&xi|(xi-stats::mean(&x))*(xi-stats::mean(&x))).sum();
    if se==0.0||x_var==0.0{return 0.0;}b/(se/crate::powers::sqrt(x_var))
}

/// Seasonal decomposition (additive): returns (trend, seasonal, residual).
pub fn decompose(data:&[f64],period:usize)->(Vec<f64>,Vec<f64>,Vec<f64>){
    let n=data.len();let trend=sma(data,period);
    // Pad trend to match data length
    let pad=(n-trend.len())/2;
    let mut full_trend=vec![trend[0];pad];full_trend.extend_from_slice(&trend);
    while full_trend.len()<n{full_trend.push(*trend.last().unwrap());}
    // Detrended
    let detrended:Vec<f64>=data.iter().zip(full_trend.iter()).map(|(&d,&t)|d-t).collect();
    // Average seasonal component
    let mut seasonal_avg=vec![0.0;period];
    for i in 0..period{let mut sum=0.0;let mut count=0;
        let mut j=i;while j<n{sum+=detrended[j];count+=1;j+=period;}seasonal_avg[i]=sum/count as f64;}
    let seasonal:Vec<f64>=(0..n).map(|i|seasonal_avg[i%period]).collect();
    let residual:Vec<f64>=(0..n).map(|i|data[i]-full_trend[i]-seasonal[i]).collect();
    (full_trend,seasonal,residual)
}

/// Forecast using simple exponential smoothing.
pub fn forecast_exp(data:&[f64],alpha:f64,horizon:usize)->Vec<f64>{
    let smoothed=ema(data,alpha);let last=*smoothed.last().unwrap();vec![last;horizon]
}

/// Ljung-Box test statistic Q.
pub fn ljung_box(data:&[f64],max_lag:usize)->f64{
    let n=data.len()as f64;let ac=acf(data,max_lag);
    n*(n+2.0)*(1..=max_lag).map(|k|ac[k]*ac[k]/(n-k as f64)).sum::<f64>()
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_sma(){let r=sma(&[1.0,2.0,3.0,4.0,5.0],3);assert_eq!(r.len(),3);assert!(approx_eq(r[0],2.0,1e-10));}
    #[test]fn test_ema(){let r=ema(&[1.0,2.0,3.0],0.5);assert!(approx_eq(r[0],1.0,1e-10));assert!(approx_eq(r[1],1.5,1e-10));}
    #[test]fn test_acf(){let r=acf(&[1.0,2.0,3.0,4.0,5.0],2);assert!(approx_eq(r[0],1.0,1e-10));}
    #[test]fn test_difference(){let r=difference(&[1.0,3.0,6.0,10.0],1);assert_eq!(r,vec![2.0,3.0,4.0]);}
}
