//! Control theory: PID, transfer functions, stability analysis.

use crate::powers;use crate::trig;use crate::constants::PI;use crate::complex::Complex;

/// PID controller state.
#[derive(Clone,Debug)]
pub struct PID{pub kp:f64,pub ki:f64,pub kd:f64,integral:f64,prev_error:f64}

impl PID{
    pub fn new(kp:f64,ki:f64,kd:f64)->Self{Self{kp,ki,kd,integral:0.0,prev_error:0.0}}
    /// Compute control output given setpoint and measured value.
    pub fn update(&mut self,setpoint:f64,measured:f64,dt:f64)->f64{
        let error=setpoint-measured;self.integral+=error*dt;
        let derivative=(error-self.prev_error)/dt;self.prev_error=error;
        self.kp*error+self.ki*self.integral+self.kd*derivative
    }
    pub fn reset(&mut self){self.integral=0.0;self.prev_error=0.0;}
}

/// Simulate PID control of first-order system: τ*dy/dt + y = K*u.
pub fn simulate_pid(pid:&mut PID,setpoint:f64,tau:f64,k:f64,dt:f64,steps:usize)->Vec<(f64,f64,f64)>{
    let mut y=0.0;let mut result=Vec::with_capacity(steps);
    for i in 0..steps{let t=i as f64*dt;let u=pid.update(setpoint,y,dt);
        y+=(k*u-y)/tau*dt;result.push((t,y,u));}result
}

/// Transfer function H(s) = num(s)/den(s) evaluated at s=jω.
pub fn transfer_function(num:&[f64],den:&[f64],omega:f64)->Complex{
    let s=Complex::new(0.0,omega);
    let eval=|coeffs:&[f64]|->Complex{let mut result=Complex::new(0.0,0.0);let mut sk=Complex::new(1.0,0.0);
        for &c in coeffs{result=Complex::new(result.re+c*sk.re,result.im+c*sk.im);
            sk=Complex::new(sk.re*s.re-sk.im*s.im,sk.re*s.im+sk.im*s.re);}result};
    let n=eval(num);let d=eval(den);
    let denom=d.re*d.re+d.im*d.im;Complex::new((n.re*d.re+n.im*d.im)/denom,(n.im*d.re-n.re*d.im)/denom)
}

/// Bode plot data: returns (frequencies, magnitudes_db, phases_deg).
pub fn bode(num:&[f64],den:&[f64],omega_start:f64,omega_end:f64,n_points:usize)->(Vec<f64>,Vec<f64>,Vec<f64>){
    let log_start=crate::exponential::log10(omega_start);let log_end=crate::exponential::log10(omega_end);
    let mut freqs=Vec::new();let mut mags=Vec::new();let mut phases=Vec::new();
    for i in 0..n_points{let log_w=log_start+(log_end-log_start)*i as f64/(n_points-1)as f64;
        let w=powers::pow(10.0,log_w);let h=transfer_function(num,den,w);
        let mag=powers::sqrt(h.re*h.re+h.im*h.im);
        freqs.push(w);mags.push(20.0*crate::exponential::log10(mag));phases.push(trig::atan2(h.im,h.re)*180.0/PI);}
    (freqs,mags,phases)
}

/// Gain margin: find frequency where phase = -180° and return gain at that frequency.
pub fn gain_margin(num:&[f64],den:&[f64])->f64{
    let(_,mags,phases)=bode(num,den,0.01,1000.0,10000);
    let mut min_phase_diff=f64::INFINITY;let mut gm_db=0.0;
    for i in 0..phases.len(){let diff=(phases[i]+180.0).abs();
        if diff<min_phase_diff{min_phase_diff=diff;gm_db=-mags[i];}}gm_db
}

/// Phase margin: find frequency where gain = 0 dB and return phase + 180°.
pub fn phase_margin(num:&[f64],den:&[f64])->f64{
    let(_,mags,phases)=bode(num,den,0.01,1000.0,10000);
    let mut min_mag=f64::INFINITY;let mut pm=0.0;
    for i in 0..mags.len(){if mags[i].abs()<min_mag{min_mag=mags[i].abs();pm=phases[i]+180.0;}}pm
}

/// Routh-Hurwitz stability criterion. Returns true if all poles in left half-plane.
pub fn routh_hurwitz(coeffs:&[f64])->bool{
    let n=coeffs.len();if n<2{return true;}
    let cols=(n+1)/2;let mut table=vec![vec![0.0;cols];n];
    for i in 0..cols{if 2*i<n{table[0][i]=coeffs[2*i];}if 2*i+1<n{table[1][i]=coeffs[2*i+1];}}
    for i in 2..n{for j in 0..cols-1{if table[i-1][0].abs()<1e-30{return false;}
        table[i][j]=(table[i-1][0]*table[i-2][j+1]-table[i-2][0]*table[i-1][j+1])/table[i-1][0];}
        if (table[i][0]<0.0)!=(table[0][0]<0.0){return false;}} // sign change
    // Check first column for sign changes
    let mut sign_changes=0;
    for i in 1..n{if table[i][0]*table[i-1][0]<0.0{sign_changes+=1;}}
    sign_changes==0
}

/// Step response of first-order system K/(τs+1).
pub fn step_response_first_order(k:f64,tau:f64,t:f64)->f64{k*(1.0-crate::exponential::exp(-t/tau))}

/// Step response of second-order system.
pub fn step_response_second_order(k:f64,wn:f64,zeta:f64,t:f64)->f64{
    if zeta>=1.0{// Overdamped
        let s1=-zeta*wn+wn*powers::sqrt(zeta*zeta-1.0);let s2=-zeta*wn-wn*powers::sqrt(zeta*zeta-1.0);
        k*(1.0+s2/(s1-s2)*crate::exponential::exp(s1*t)-s1/(s1-s2)*crate::exponential::exp(s2*t))
    }else{// Underdamped
        let wd=wn*powers::sqrt(1.0-zeta*zeta);
        k*(1.0-crate::exponential::exp(-zeta*wn*t)/powers::sqrt(1.0-zeta*zeta)*trig::sin(wd*t+trig::acos(zeta)))
    }
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_pid(){let mut pid=PID::new(1.0,0.1,0.01);let u=pid.update(10.0,0.0,0.1);assert!(u>0.0);}
    #[test]fn test_step_first(){assert!(approx_eq(step_response_first_order(1.0,1.0,0.0),0.0,1e-10));}
    #[test]fn test_routh(){assert!(routh_hurwitz(&[1.0,3.0,3.0,1.0]));} // s³+3s²+3s+1 = (s+1)³
}
