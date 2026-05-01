//! Random sampling from distributions.

use crate::powers;
use crate::exponential;
use crate::trig;
use crate::constants;

/// Linear congruential generator state.
pub struct Rng{state:u64}
impl Rng{
    pub fn new(seed:u64)->Self{Self{state:seed}}
    pub fn next_u64(&mut self)->u64{self.state=self.state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);self.state}
    /// Uniform f64 in [0, 1).
    pub fn uniform(&mut self)->f64{(self.next_u64()>>11)as f64/9007199254740992.0}
    /// Uniform f64 in [a, b).
    pub fn uniform_range(&mut self,a:f64,b:f64)->f64{a+self.uniform()*(b-a)}
    /// Uniform integer in [a, b].
    pub fn int_range(&mut self,a:i64,b:i64)->i64{a+(self.next_u64()%((b-a+1)as u64))as i64}
}

/// Box-Muller transform: generate standard normal sample.
pub fn box_muller(rng:&mut Rng)->(f64,f64){
    let u1=rng.uniform().max(1e-30);let u2=rng.uniform();
    let r=powers::sqrt(-2.0*exponential::ln(u1));let theta=constants::TAU*u2;
    (r*trig::cos(theta),r*trig::sin(theta))
}

/// Sample from normal distribution N(mean, sigma).
pub fn sample_normal(rng:&mut Rng,mean:f64,sigma:f64)->f64{let(z,_)=box_muller(rng);mean+sigma*z}

/// Sample n values from normal distribution.
pub fn sample_normal_vec(rng:&mut Rng,mean:f64,sigma:f64,n:usize)->Vec<f64>{
    (0..n).map(|_|sample_normal(rng,mean,sigma)).collect()
}

/// Sample from exponential distribution with rate λ.
pub fn sample_exponential(rng:&mut Rng,lambda:f64)->f64{-exponential::ln(rng.uniform().max(1e-30))/lambda}

/// Sample from Poisson distribution with mean λ (Knuth's algorithm).
pub fn sample_poisson(rng:&mut Rng,lambda:f64)->u64{
    let l=exponential::exp(-lambda);let mut k=0u64;let mut p=1.0;
    loop{k+=1;p*=rng.uniform();if p<l{return k-1;}}
}

/// Sample from geometric distribution with probability p.
pub fn sample_geometric(rng:&mut Rng,p:f64)->u64{
    (exponential::ln(rng.uniform().max(1e-30))/exponential::ln(1.0-p)).ceil()as u64
}

/// Sample from Bernoulli distribution.
pub fn sample_bernoulli(rng:&mut Rng,p:f64)->bool{rng.uniform()<p}

/// Sample from binomial distribution B(n, p).
pub fn sample_binomial(rng:&mut Rng,n:u64,p:f64)->u64{(0..n).filter(|_|sample_bernoulli(rng,p)).count()as u64}

/// Sample from uniform discrete [a, b].
pub fn sample_discrete_uniform(rng:&mut Rng,a:i64,b:i64)->i64{rng.int_range(a,b)}

/// Sample from categorical distribution given weights.
pub fn sample_categorical(rng:&mut Rng,weights:&[f64])->usize{
    let total:f64=weights.iter().sum();let mut r=rng.uniform()*total;
    for(i,&w) in weights.iter().enumerate(){r-=w;if r<=0.0{return i;}}weights.len()-1
}

/// Reservoir sampling: select k items uniformly from a stream.
pub fn reservoir_sample(rng:&mut Rng,data:&[f64],k:usize)->Vec<f64>{
    let mut reservoir=data[..k].to_vec();
    for i in k..data.len(){let j=rng.int_range(0,i as i64)as usize;if j<k{reservoir[j]=data[i];}}reservoir
}

/// Fisher-Yates shuffle.
pub fn shuffle<T:Clone>(rng:&mut Rng,data:&mut [T]){
    let n=data.len();for i in(1..n).rev(){let j=rng.int_range(0,i as i64)as usize;data.swap(i,j);}
}

/// Rejection sampling from target density using uniform proposal.
pub fn rejection_sample(rng:&mut Rng,target:&dyn Fn(f64)->f64,x_min:f64,x_max:f64,y_max:f64,n:usize)->Vec<f64>{
    let mut samples=Vec::with_capacity(n);
    while samples.len()<n{let x=rng.uniform_range(x_min,x_max);let y=rng.uniform_range(0.0,y_max);
        if y<=target(x){samples.push(x);}}samples
}

/// Monte Carlo estimate of E[f(X)] where X ~ Uniform[a,b].
pub fn monte_carlo_expectation(rng:&mut Rng,f:&dyn Fn(f64)->f64,a:f64,b:f64,n:usize)->f64{
    let mut sum=0.0;for _ in 0..n{sum+=f(rng.uniform_range(a,b));}sum/n as f64
}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_normal(){let mut rng=Rng::new(42);let samples=sample_normal_vec(&mut rng,0.0,1.0,10000);
        let mean:f64=samples.iter().sum::<f64>()/10000.0;assert!(mean.abs()<0.1);}
    #[test]fn test_uniform(){let mut rng=Rng::new(42);let samples:Vec<f64>=(0..1000).map(|_|rng.uniform()).collect();
        assert!(samples.iter().all(|&x|x>=0.0&&x<1.0));}
    #[test]fn test_poisson(){let mut rng=Rng::new(42);let samples:Vec<u64>=(0..10000).map(|_|sample_poisson(&mut rng,5.0)).collect();
        let mean=samples.iter().sum::<u64>()as f64/10000.0;assert!((mean-5.0).abs()<0.3);}
    #[test]fn test_shuffle(){let mut rng=Rng::new(42);let mut v=vec![1,2,3,4,5];shuffle(&mut rng,&mut v);
        v.sort();assert_eq!(v,vec![1,2,3,4,5]);}// all elements preserved
}
