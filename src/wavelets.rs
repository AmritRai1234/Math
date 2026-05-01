//! Wavelet transforms: Haar, Daubechies.

/// Haar wavelet transform (in-place, signal length must be power of 2).
pub fn haar_forward(data:&[f64])->Vec<f64>{
    let mut result=data.to_vec();let mut n=result.len();
    while n>1{let mut temp=vec![0.0;n];
        for i in 0..n/2{temp[i]=(result[2*i]+result[2*i+1])/crate::constants::SQRT2;
            temp[n/2+i]=(result[2*i]-result[2*i+1])/crate::constants::SQRT2;}
        for i in 0..n{result[i]=temp[i];}n/=2;}
    result
}

/// Inverse Haar wavelet transform.
pub fn haar_inverse(data:&[f64])->Vec<f64>{
    let mut result=data.to_vec();let len=result.len();let mut n=2;
    while n<=len{let mut temp=vec![0.0;n];
        for i in 0..n/2{temp[2*i]=(result[i]+result[n/2+i])/crate::constants::SQRT2;
            temp[2*i+1]=(result[i]-result[n/2+i])/crate::constants::SQRT2;}
        for i in 0..n{result[i]=temp[i];}n*=2;}
    result
}

/// Daubechies D4 wavelet forward transform.
pub fn daubechies4_forward(data:&[f64])->Vec<f64>{
    let c0=(1.0+crate::powers::sqrt(3.0))/(4.0*crate::constants::SQRT2);
    let c1=(3.0+crate::powers::sqrt(3.0))/(4.0*crate::constants::SQRT2);
    let c2=(3.0-crate::powers::sqrt(3.0))/(4.0*crate::constants::SQRT2);
    let c3=(1.0-crate::powers::sqrt(3.0))/(4.0*crate::constants::SQRT2);
    let mut result=data.to_vec();let mut n=result.len();
    while n>=4{let mut temp=vec![0.0;n];
        for i in 0..n/2{let j=2*i;
            temp[i]=c0*result[j%n]+c1*result[(j+1)%n]+c2*result[(j+2)%n]+c3*result[(j+3)%n];
            temp[n/2+i]=c3*result[j%n]-c2*result[(j+1)%n]+c1*result[(j+2)%n]-c0*result[(j+3)%n];}
        for i in 0..n{result[i]=temp[i];}n/=2;}
    result
}

/// Wavelet thresholding (hard).
pub fn hard_threshold(coeffs:&mut [f64],threshold:f64){
    for c in coeffs.iter_mut(){if crate::arithmetic::abs(*c)<threshold{*c=0.0;}}
}

/// Wavelet thresholding (soft).
pub fn soft_threshold(coeffs:&mut [f64],threshold:f64){
    for c in coeffs.iter_mut(){if crate::arithmetic::abs(*c)<threshold{*c=0.0;}
        else if *c>0.0{*c-=threshold;}else{*c+=threshold;}}
}

/// Wavelet denoising using Haar.
pub fn denoise_haar(data:&[f64],threshold:f64)->Vec<f64>{
    let mut coeffs=haar_forward(data);
    soft_threshold(&mut coeffs[1..],threshold); // keep scaling coeff
    haar_inverse(&coeffs)
}

/// Multi-resolution analysis: returns approximation and detail coefficients at each level.
pub fn mra(data:&[f64],levels:usize)->Vec<(Vec<f64>,Vec<f64>)>{
    let mut result=Vec::new();let mut current=data.to_vec();
    for _ in 0..levels{if current.len()<2{break;}
        let n=current.len();let mut approx=Vec::new();let mut detail=Vec::new();
        for i in 0..n/2{approx.push((current[2*i]+current[2*i+1])/crate::constants::SQRT2);
            detail.push((current[2*i]-current[2*i+1])/crate::constants::SQRT2);}
        result.push((approx.clone(),detail));current=approx;}
    result
}

/// Energy of wavelet coefficients.
pub fn wavelet_energy(coeffs:&[f64])->f64{coeffs.iter().map(|&c|c*c).sum()}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_haar_roundtrip(){let data=[1.0,2.0,3.0,4.0];let fwd=haar_forward(&data);let inv=haar_inverse(&fwd);
        for i in 0..4{assert!(approx_eq(inv[i],data[i],1e-10));}}
    #[test]fn test_haar_energy(){let data=[1.0,2.0,3.0,4.0];let coeffs=haar_forward(&data);
        let e_orig:f64=data.iter().map(|&x|x*x).sum();let e_coeffs:f64=coeffs.iter().map(|&x|x*x).sum();
        assert!(approx_eq(e_orig,e_coeffs,1e-10));}// Parseval's theorem
    #[test]fn test_mra(){let data=[1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0];let levels=mra(&data,3);assert_eq!(levels.len(),3);}
}
