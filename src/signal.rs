//! Signal processing basics: convolution, correlation, DFT.

use crate::trig;use crate::constants::PI;

/// Discrete Fourier Transform.
pub fn dft(signal:&[f64])->Vec<(f64,f64)>{
    let n=signal.len();
    (0..n).map(|k|{let mut re=0.0;let mut im=0.0;
        for t in 0..n{let angle=-2.0*PI*k as f64*t as f64/n as f64;
            re+=signal[t]*trig::cos(angle);im+=signal[t]*trig::sin(angle);}(re,im)}).collect()
}

/// Inverse DFT.
pub fn idft(spectrum:&[(f64,f64)])->Vec<f64>{
    let n=spectrum.len();
    (0..n).map(|t|{let mut re=0.0;
        for k in 0..n{let angle=2.0*PI*k as f64*t as f64/n as f64;
            re+=spectrum[k].0*trig::cos(angle)-spectrum[k].1*trig::sin(angle);}re/n as f64}).collect()
}

/// 1D convolution.
pub fn convolve(a:&[f64],b:&[f64])->Vec<f64>{
    let n=a.len()+b.len()-1;
    (0..n).map(|i|{let mut sum=0.0;for j in 0..b.len(){if i>=j&&i-j<a.len(){sum+=a[i-j]*b[j];}}sum}).collect()
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_dft_dc(){let s=vec![1.0,1.0,1.0,1.0];let f=dft(&s);assert!(approx_eq(f[0].0,4.0,1e-10));}
    #[test]fn test_convolve(){let c=convolve(&[1.0,2.0],&[3.0,4.0]);assert_eq!(c.len(),3);assert!(approx_eq(c[0],3.0,1e-10));}
}
