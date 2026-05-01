//! Color math: color space conversions, perceptual distance.

use crate::powers;use crate::trig;use crate::constants::PI;

/// RGB (0-255) to HSL (h:0-360, s:0-1, l:0-1).
pub fn rgb_to_hsl(r:u8,g:u8,b:u8)->(f64,f64,f64){
    let(r,g,b)=(r as f64/255.0,g as f64/255.0,b as f64/255.0);
    let max=r.max(g).max(b);let min=r.min(g).min(b);let l=(max+min)/2.0;
    if(max-min).abs()<1e-10{return(0.0,0.0,l);}
    let d=max-min;let s=if l>0.5{d/(2.0-max-min)}else{d/(max+min)};
    let h=if(max-r).abs()<1e-10{(g-b)/d+(if g<b{6.0}else{0.0})}
        else if(max-g).abs()<1e-10{(b-r)/d+2.0}else{(r-g)/d+4.0};
    (h*60.0,s,l)
}

/// HSL to RGB.
pub fn hsl_to_rgb(h:f64,s:f64,l:f64)->(u8,u8,u8){
    if s==0.0{let v=(l*255.0)as u8;return(v,v,v);}
    let hue2rgb=|p:f64,q:f64,mut t:f64|->f64{if t<0.0{t+=1.0;}if t>1.0{t-=1.0;}
        if t<1.0/6.0{p+(q-p)*6.0*t}else if t<0.5{q}else if t<2.0/3.0{p+(q-p)*(2.0/3.0-t)*6.0}else{p}};
    let q=if l<0.5{l*(1.0+s)}else{l+s-l*s};let p=2.0*l-q;let h=h/360.0;
    ((hue2rgb(p,q,h+1.0/3.0)*255.0)as u8,(hue2rgb(p,q,h)*255.0)as u8,(hue2rgb(p,q,h-1.0/3.0)*255.0)as u8)
}

/// RGB to HSV (h:0-360, s:0-1, v:0-1).
pub fn rgb_to_hsv(r:u8,g:u8,b:u8)->(f64,f64,f64){
    let(r,g,b)=(r as f64/255.0,g as f64/255.0,b as f64/255.0);
    let max=r.max(g).max(b);let min=r.min(g).min(b);let d=max-min;
    let v=max;let s=if max==0.0{0.0}else{d/max};
    if d<1e-10{return(0.0,s,v);}
    let h=if(max-r).abs()<1e-10{(g-b)/d+(if g<b{6.0}else{0.0})}
        else if(max-g).abs()<1e-10{(b-r)/d+2.0}else{(r-g)/d+4.0};
    (h*60.0,s,v)
}

/// RGB to CIE XYZ.
pub fn rgb_to_xyz(r:u8,g:u8,b:u8)->(f64,f64,f64){
    let linearize=|c:f64|if c>0.04045{powers::pow((c+0.055)/1.055,2.4)}else{c/12.92};
    let(r,g,b)=(linearize(r as f64/255.0),linearize(g as f64/255.0),linearize(b as f64/255.0));
    (0.4124564*r+0.3575761*g+0.1804375*b,0.2126729*r+0.7151522*g+0.0721750*b,0.0193339*r+0.1191920*g+0.9503041*b)
}

/// CIE XYZ to CIELAB.
pub fn xyz_to_lab(x:f64,y:f64,z:f64)->(f64,f64,f64){
    let(xn,yn,zn)=(0.95047,1.0,1.08883); // D65 illuminant
    let f=|t:f64|if t>0.008856{powers::pow(t,1.0/3.0)}else{7.787*t+16.0/116.0};
    let l=116.0*f(y/yn)-16.0;let a=500.0*(f(x/xn)-f(y/yn));let b=200.0*(f(y/yn)-f(z/zn));
    (l,a,b)
}

/// RGB to CIELAB.
pub fn rgb_to_lab(r:u8,g:u8,b:u8)->(f64,f64,f64){let(x,y,z)=rgb_to_xyz(r,g,b);xyz_to_lab(x,y,z)}

/// CIEDE2000 color difference.
pub fn ciede2000(l1:f64,a1:f64,b1:f64,l2:f64,a2:f64,b2:f64)->f64{
    let dl=l2-l1;let ml=(l1+l2)/2.0;
    let c1=powers::sqrt(a1*a1+b1*b1);let c2=powers::sqrt(a2*a2+b2*b2);let mc=(c1+c2)/2.0;
    let mc7=powers::pow(mc,7.0);let g=0.5*(1.0-powers::sqrt(mc7/(mc7+powers::pow(25.0,7.0))));
    let a1p=a1*(1.0+g);let a2p=a2*(1.0+g);
    let c1p=powers::sqrt(a1p*a1p+b1*b1);let c2p=powers::sqrt(a2p*a2p+b2*b2);
    let dc=c2p-c1p;let mcp=(c1p+c2p)/2.0;
    let h1=trig::atan2(b1,a1p)*180.0/PI;let h1=if h1<0.0{h1+360.0}else{h1};
    let h2=trig::atan2(b2,a2p)*180.0/PI;let h2=if h2<0.0{h2+360.0}else{h2};
    let dh=if(h1-h2).abs()<=180.0{h2-h1}else if h2<=h1{h2-h1+360.0}else{h2-h1-360.0};
    let dhp=2.0*powers::sqrt(c1p*c2p)*trig::sin(dh*PI/360.0);
    let mhp=if(h1-h2).abs()<=180.0{(h1+h2)/2.0}else{(h1+h2+360.0)/2.0};
    let t=1.0-0.17*trig::cos((mhp-30.0)*PI/180.0)+0.24*trig::cos((2.0*mhp)*PI/180.0)
        +0.32*trig::cos((3.0*mhp+6.0)*PI/180.0)-0.20*trig::cos((4.0*mhp-63.0)*PI/180.0);
    let sl=1.0+0.015*(ml-50.0)*(ml-50.0)/powers::sqrt(20.0+(ml-50.0)*(ml-50.0));
    let sc=1.0+0.045*mcp;let sh=1.0+0.015*mcp*t;
    let mcp7=powers::pow(mcp,7.0);
    let rt=-2.0*powers::sqrt(mcp7/(mcp7+powers::pow(25.0,7.0)))*trig::sin(60.0*crate::exponential::exp(-((mhp-275.0)/25.0)*((mhp-275.0)/25.0))*PI/180.0);
    powers::sqrt((dl/sl)*(dl/sl)+(dc/sc)*(dc/sc)+(dhp/sh)*(dhp/sh)+rt*(dc/sc)*(dhp/sh))
}

/// Hex string to RGB.
pub fn hex_to_rgb(hex:&str)->(u8,u8,u8){
    let hex=hex.trim_start_matches('#');
    let r=u8::from_str_radix(&hex[0..2],16).unwrap_or(0);
    let g=u8::from_str_radix(&hex[2..4],16).unwrap_or(0);
    let b=u8::from_str_radix(&hex[4..6],16).unwrap_or(0);(r,g,b)
}

/// RGB to hex string.
pub fn rgb_to_hex(r:u8,g:u8,b:u8)->String{format!("#{:02X}{:02X}{:02X}",r,g,b)}

/// Linear interpolation between two colors.
pub fn lerp_color(r1:u8,g1:u8,b1:u8,r2:u8,g2:u8,b2:u8,t:f64)->(u8,u8,u8){
    let lerp=|a:u8,b:u8|(a as f64+(b as f64-a as f64)*t)as u8;
    (lerp(r1,r2),lerp(g1,g2),lerp(b1,b2))
}

/// Relative luminance (WCAG).
pub fn relative_luminance(r:u8,g:u8,b:u8)->f64{
    let lin=|c:u8|{let c=c as f64/255.0;if c<=0.03928{c/12.92}else{powers::pow((c+0.055)/1.055,2.4)}};
    0.2126*lin(r)+0.7152*lin(g)+0.0722*lin(b)
}

/// WCAG contrast ratio.
pub fn contrast_ratio(r1:u8,g1:u8,b1:u8,r2:u8,g2:u8,b2:u8)->f64{
    let l1=relative_luminance(r1,g1,b1);let l2=relative_luminance(r2,g2,b2);
    let(lighter,darker)=if l1>l2{(l1,l2)}else{(l2,l1)};(lighter+0.05)/(darker+0.05)
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_hsl_roundtrip(){let(h,s,l)=rgb_to_hsl(128,64,192);let(r,g,b)=hsl_to_rgb(h,s,l);
        assert!((r as i16-128).abs()<=1&&(g as i16-64).abs()<=1&&(b as i16-192).abs()<=1);}
    #[test]fn test_hex(){assert_eq!(hex_to_rgb("#FF8000"),(255,128,0));assert_eq!(rgb_to_hex(255,128,0),"#FF8000");}
    #[test]fn test_contrast(){let c=contrast_ratio(255,255,255,0,0,0);assert!(c>20.0);} // black/white
    #[test]fn test_ciede_same(){assert!(approx_eq(ciede2000(50.0,0.0,0.0,50.0,0.0,0.0),0.0,1e-6));}
}
