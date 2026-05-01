//! Image processing kernels: convolution, edge detection, morphological operations.

/// 2D convolution of image with kernel.
pub fn convolve2d(image:&[Vec<f64>],kernel:&[Vec<f64>])->Vec<Vec<f64>>{
    let(ih,iw)=(image.len(),image[0].len());let(kh,kw)=(kernel.len(),kernel[0].len());
    let(ph,pw)=(kh/2,kw/2);
    (0..ih).map(|i|(0..iw).map(|j|{let mut sum=0.0;
        for ki in 0..kh{for kj in 0..kw{let ii=i as i32+ki as i32-ph as i32;let jj=j as i32+kj as i32-pw as i32;
            if ii>=0&&ii<ih as i32&&jj>=0&&jj<iw as i32{sum+=image[ii as usize][jj as usize]*kernel[ki][kj];}}}sum}).collect()).collect()
}

/// Sobel edge detection (returns gradient magnitude).
pub fn sobel(image:&[Vec<f64>])->Vec<Vec<f64>>{
    let gx=convolve2d(image,&[vec![-1.0,0.0,1.0],vec![-2.0,0.0,2.0],vec![-1.0,0.0,1.0]]);
    let gy=convolve2d(image,&[vec![-1.0,-2.0,-1.0],vec![0.0,0.0,0.0],vec![1.0,2.0,1.0]]);
    gx.iter().zip(gy.iter()).map(|(rx,ry)|rx.iter().zip(ry.iter()).map(|(&x,&y)|crate::powers::sqrt(x*x+y*y)).collect()).collect()
}

/// Gaussian blur kernel.
pub fn gaussian_kernel(size:usize,sigma:f64)->Vec<Vec<f64>>{
    let center=size as f64/2.0;let mut kernel=vec![vec![0.0;size];size];let mut sum=0.0;
    for i in 0..size{for j in 0..size{let x=i as f64-center+0.5;let y=j as f64-center+0.5;
        kernel[i][j]=crate::exponential::exp(-(x*x+y*y)/(2.0*sigma*sigma));sum+=kernel[i][j];}}
    for row in &mut kernel{for v in row{*v/=sum;}}kernel
}

/// Gaussian blur.
pub fn gaussian_blur(image:&[Vec<f64>],size:usize,sigma:f64)->Vec<Vec<f64>>{
    convolve2d(image,&gaussian_kernel(size,sigma))
}

/// Sharpen kernel.
pub fn sharpen(image:&[Vec<f64>])->Vec<Vec<f64>>{
    convolve2d(image,&[vec![0.0,-1.0,0.0],vec![-1.0,5.0,-1.0],vec![0.0,-1.0,0.0]])
}

/// Laplacian edge detection.
pub fn laplacian(image:&[Vec<f64>])->Vec<Vec<f64>>{
    convolve2d(image,&[vec![0.0,1.0,0.0],vec![1.0,-4.0,1.0],vec![0.0,1.0,0.0]])
}

/// Prewitt edge detection.
pub fn prewitt(image:&[Vec<f64>])->Vec<Vec<f64>>{
    let gx=convolve2d(image,&[vec![-1.0,0.0,1.0],vec![-1.0,0.0,1.0],vec![-1.0,0.0,1.0]]);
    let gy=convolve2d(image,&[vec![-1.0,-1.0,-1.0],vec![0.0,0.0,0.0],vec![1.0,1.0,1.0]]);
    gx.iter().zip(gy.iter()).map(|(rx,ry)|rx.iter().zip(ry.iter()).map(|(&x,&y)|crate::powers::sqrt(x*x+y*y)).collect()).collect()
}

/// Morphological erosion (binary image, structuring element is ones).
pub fn erode(image:&[Vec<u8>],se_size:usize)->Vec<Vec<u8>>{
    let(h,w)=(image.len(),image[0].len());let pad=se_size/2;
    (0..h).map(|i|(0..w).map(|j|{for di in 0..se_size{for dj in 0..se_size{
        let ni=i as i32+di as i32-pad as i32;let nj=j as i32+dj as i32-pad as i32;
        if ni<0||ni>=h as i32||nj<0||nj>=w as i32||image[ni as usize][nj as usize]==0{return 0;}}}1}).collect()).collect()
}

/// Morphological dilation (binary image).
pub fn dilate(image:&[Vec<u8>],se_size:usize)->Vec<Vec<u8>>{
    let(h,w)=(image.len(),image[0].len());let pad=se_size/2;
    (0..h).map(|i|(0..w).map(|j|{for di in 0..se_size{for dj in 0..se_size{
        let ni=i as i32+di as i32-pad as i32;let nj=j as i32+dj as i32-pad as i32;
        if ni>=0&&ni<h as i32&&nj>=0&&nj<w as i32&&image[ni as usize][nj as usize]==1{return 1;}}}0}).collect()).collect()
}

/// Morphological opening (erosion then dilation).
pub fn opening(image:&[Vec<u8>],se_size:usize)->Vec<Vec<u8>>{dilate(&erode(image,se_size),se_size)}
/// Morphological closing (dilation then erosion).
pub fn closing(image:&[Vec<u8>],se_size:usize)->Vec<Vec<u8>>{erode(&dilate(image,se_size),se_size)}

/// Image histogram (256 bins for grayscale).
pub fn histogram(image:&[Vec<u8>])->Vec<u32>{
    let mut hist=vec![0u32;256];for row in image{for &p in row{hist[p as usize]+=1;}}hist
}

/// Otsu's threshold.
pub fn otsu_threshold(hist:&[u32])->u8{
    let total:u32=hist.iter().sum();let mut sum_total:f64=hist.iter().enumerate().map(|(i,&c)|i as f64*c as f64).sum();
    let mut sum_bg=0.0;let mut w_bg=0u32;let mut max_var=0.0;let mut threshold=0u8;
    for t in 0..256{w_bg+=hist[t];if w_bg==0{continue;}let w_fg=total-w_bg;if w_fg==0{break;}
        sum_bg+=t as f64*hist[t]as f64;let mean_bg=sum_bg/w_bg as f64;let mean_fg=(sum_total-sum_bg)/w_fg as f64;
        let var=w_bg as f64*w_fg as f64*(mean_bg-mean_fg)*(mean_bg-mean_fg);
        if var>max_var{max_var=var;threshold=t as u8;}}threshold
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_convolve_identity(){let img=vec![vec![1.0,2.0,3.0],vec![4.0,5.0,6.0],vec![7.0,8.0,9.0]];
        let kernel=vec![vec![0.0,0.0,0.0],vec![0.0,1.0,0.0],vec![0.0,0.0,0.0]];
        let result=convolve2d(&img,&kernel);assert!(approx_eq(result[1][1],5.0,1e-10));}
    #[test]fn test_sobel(){let img=vec![vec![0.0,0.0,0.0],vec![0.0,1.0,0.0],vec![0.0,0.0,0.0]];
        let edges=sobel(&img);assert!(edges[1][1]>=0.0);}
    #[test]fn test_dilate(){let img=vec![vec![0,0,0],vec![0,1,0],vec![0,0,0]];
        let d=dilate(&img,3);assert_eq!(d[0][0],1);}
}
