//! Signal processing: windowing, filters, spectral analysis, correlation.

use crate::trig;use crate::constants::PI;use crate::powers;

// ─── Window Functions ───

/// Hann window.
pub fn hann(n:usize)->Vec<f64>{(0..n).map(|i|0.5*(1.0-trig::cos(2.0*PI*i as f64/(n-1)as f64))).collect()}
/// Hamming window.
pub fn hamming(n:usize)->Vec<f64>{(0..n).map(|i|0.54-0.46*trig::cos(2.0*PI*i as f64/(n-1)as f64)).collect()}
/// Blackman window.
pub fn blackman(n:usize)->Vec<f64>{(0..n).map(|i|0.42-0.5*trig::cos(2.0*PI*i as f64/(n-1)as f64)+0.08*trig::cos(4.0*PI*i as f64/(n-1)as f64)).collect()}
/// Bartlett (triangular) window.
pub fn bartlett(n:usize)->Vec<f64>{(0..n).map(|i|1.0-crate::arithmetic::abs(2.0*i as f64/(n-1)as f64-1.0)).collect()}
/// Kaiser window (approximate).
pub fn kaiser(n:usize,beta:f64)->Vec<f64>{
    let i0b=bessel_i0(beta);
    (0..n).map(|i|{let x=2.0*i as f64/(n-1)as f64-1.0;bessel_i0(beta*powers::sqrt(1.0-x*x))/i0b}).collect()
}
fn bessel_i0(x:f64)->f64{let mut sum=1.0;let mut term=1.0;for k in 1..20{term*=(x/2.0)*(x/2.0)/(k*k)as f64;sum+=term;}sum}

/// Apply window to signal.
pub fn apply_window(signal:&[f64],window:&[f64])->Vec<f64>{
    signal.iter().zip(window.iter()).map(|(&s,&w)|s*w).collect()
}

// ─── Filters ───

/// FIR filter (finite impulse response): y[n] = Σ h[k] * x[n-k].
pub fn fir_filter(signal:&[f64],coeffs:&[f64])->Vec<f64>{
    let n=signal.len();let m=coeffs.len();
    (0..n).map(|i|{let mut sum=0.0;for k in 0..m{if i>=k{sum+=coeffs[k]*signal[i-k];}}sum}).collect()
}

/// IIR filter: y[n] = Σ b[k]*x[n-k] - Σ a[k]*y[n-k].
pub fn iir_filter(signal:&[f64],b:&[f64],a:&[f64])->Vec<f64>{
    let n=signal.len();let mut y=vec![0.0;n];
    for i in 0..n{let mut sum=0.0;
        for k in 0..b.len(){if i>=k{sum+=b[k]*signal[i-k];}}
        for k in 1..a.len(){if i>=k{sum-=a[k]*y[i-k];}}
        y[i]=sum/a[0];}y
}

/// Moving average filter.
pub fn moving_average_filter(signal:&[f64],window:usize)->Vec<f64>{
    let coeffs=vec![1.0/window as f64;window];fir_filter(signal,&coeffs)
}

/// Simple low-pass RC filter (first-order IIR).
pub fn lowpass_rc(signal:&[f64],alpha:f64)->Vec<f64>{
    let b=[alpha];let a=[1.0,-(1.0-alpha)];iir_filter(signal,&b,&a)
}

/// Simple high-pass filter.
pub fn highpass(signal:&[f64],alpha:f64)->Vec<f64>{
    let lp=lowpass_rc(signal,alpha);signal.iter().zip(lp.iter()).map(|(&s,&l)|s-l).collect()
}

// ─── Spectral Analysis ───

/// Power spectral density (magnitude squared of DFT).
pub fn power_spectrum(signal:&[f64])->Vec<f64>{
    let n=signal.len();
    (0..n).map(|k|{let mut re=0.0;let mut im=0.0;
        for t in 0..n{let angle=-2.0*PI*k as f64*t as f64/n as f64;
            re+=signal[t]*trig::cos(angle);im+=signal[t]*trig::sin(angle);}re*re+im*im}).collect()
}

/// Cross-correlation of two signals.
pub fn cross_correlation(a:&[f64],b:&[f64])->Vec<f64>{
    let n=a.len().max(b.len());let mut result=vec![0.0;2*n-1];
    for lag in 0..2*n-1{let offset=lag as i64-n as i64+1;
        for i in 0..a.len(){let j=i as i64-offset;
            if j>=0&&(j as usize)<b.len(){result[lag]+=a[i]*b[j as usize];}}}result
}

/// Autocorrelation (cross-correlation with self).
pub fn autocorrelation(signal:&[f64])->Vec<f64>{cross_correlation(signal,signal)}

/// Spectral centroid.
pub fn spectral_centroid(magnitudes:&[f64],sample_rate:f64)->f64{
    let n=magnitudes.len();
    let num:f64=(0..n).map(|i|i as f64*sample_rate/n as f64*magnitudes[i]).sum();
    let den:f64=magnitudes.iter().sum();if den==0.0{0.0}else{num/den}
}

/// Zero-crossing rate.
pub fn zero_crossing_rate(signal:&[f64])->f64{
    let crossings:usize=signal.windows(2).filter(|w|(w[0]>=0.0)!=(w[1]>=0.0)).count();
    crossings as f64/(signal.len()-1)as f64
}

/// RMS (root mean square) energy.
pub fn rms(signal:&[f64])->f64{powers::sqrt(signal.iter().map(|&x|x*x).sum::<f64>()/signal.len()as f64)}

/// Signal-to-noise ratio in dB.
pub fn snr_db(signal:&[f64],noise:&[f64])->f64{
    let sig_power:f64=signal.iter().map(|&x|x*x).sum::<f64>()/signal.len()as f64;
    let noise_power:f64=noise.iter().map(|&x|x*x).sum::<f64>()/noise.len()as f64;
    10.0*crate::exponential::log10(sig_power/noise_power)
}

/// Peak signal-to-noise ratio.
pub fn psnr(original:&[f64],noisy:&[f64],max_val:f64)->f64{
    let mse:f64=original.iter().zip(noisy.iter()).map(|(&a,&b)|(a-b)*(a-b)).sum::<f64>()/original.len()as f64;
    10.0*crate::exponential::log10(max_val*max_val/mse)
}

/// Downsample signal by factor.
pub fn downsample(signal:&[f64],factor:usize)->Vec<f64>{signal.iter().step_by(factor).cloned().collect()}

/// Upsample signal by factor (zero-insertion).
pub fn upsample(signal:&[f64],factor:usize)->Vec<f64>{
    let mut result=vec![0.0;signal.len()*factor];
    for(i,&s) in signal.iter().enumerate(){result[i*factor]=s;}result
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_hann(){let w=hann(5);assert!(approx_eq(w[0],0.0,1e-10));assert!(approx_eq(w[2],1.0,1e-10));}
    #[test]fn test_fir(){let signal=[1.0,2.0,3.0,4.0,5.0];let coeffs=[0.5,0.5];
        let filtered=fir_filter(&signal,&coeffs);assert!(approx_eq(filtered[1],1.5,1e-10));}
    #[test]fn test_zcr(){let signal=[1.0,-1.0,1.0,-1.0];assert!(approx_eq(zero_crossing_rate(&signal),1.0,1e-10));}
    #[test]fn test_rms(){assert!(approx_eq(rms(&[3.0,4.0]),powers::sqrt(12.5),1e-10));}
}
