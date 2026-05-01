//! Information theory: entropy, KL divergence, mutual information.

use crate::exponential;
use crate::arithmetic::abs;

/// Shannon entropy H(X) = -Σ p(x) log₂(p(x)).
pub fn entropy(probs:&[f64])->f64{
    probs.iter().filter(|&&p|p>0.0).map(|&p|-p*exponential::log2(p)).sum()
}

/// Binary entropy H(p).
pub fn binary_entropy(p:f64)->f64{
    if p<=0.0||p>=1.0{return 0.0;}
    -p*exponential::log2(p)-(1.0-p)*exponential::log2(1.0-p)
}

/// Joint entropy H(X,Y) from joint probability table.
pub fn joint_entropy(joint:&[Vec<f64>])->f64{
    joint.iter().flat_map(|row|row.iter()).filter(|&&p|p>0.0).map(|&p|-p*exponential::log2(p)).sum()
}

/// Conditional entropy H(X|Y) = H(X,Y) - H(Y).
pub fn conditional_entropy(joint:&[Vec<f64>])->f64{
    let n=joint[0].len();let py:Vec<f64>=(0..n).map(|j|joint.iter().map(|row|row[j]).sum()).collect();
    joint_entropy(joint)-entropy(&py)
}

/// KL divergence D_KL(P || Q) = Σ p(x) log(p(x)/q(x)).
pub fn kl_divergence(p:&[f64],q:&[f64])->f64{
    p.iter().zip(q.iter()).filter(|(&pi,&qi)|pi>0.0&&qi>0.0)
        .map(|(&pi,&qi)|pi*exponential::ln(pi/qi)).sum()
}

/// Jensen-Shannon divergence.
pub fn js_divergence(p:&[f64],q:&[f64])->f64{
    let m:Vec<f64>=p.iter().zip(q.iter()).map(|(&pi,&qi)|(pi+qi)/2.0).collect();
    (kl_divergence(p,&m)+kl_divergence(q,&m))/2.0
}

/// Cross-entropy H(P, Q) = -Σ p(x) log(q(x)).
pub fn cross_entropy(p:&[f64],q:&[f64])->f64{
    p.iter().zip(q.iter()).filter(|(&_pi,&qi)|qi>0.0)
        .map(|(&pi,&qi)|-pi*exponential::ln(qi)).sum()
}

/// Mutual information I(X;Y) = H(X) + H(Y) - H(X,Y).
pub fn mutual_information(joint:&[Vec<f64>])->f64{
    let m=joint.len();let n=joint[0].len();
    let px:Vec<f64>=(0..m).map(|i|joint[i].iter().sum()).collect();
    let py:Vec<f64>=(0..n).map(|j|joint.iter().map(|row|row[j]).sum()).collect();
    entropy(&px)+entropy(&py)-joint_entropy(joint)
}

/// Normalized mutual information NMI.
pub fn normalized_mutual_information(joint:&[Vec<f64>])->f64{
    let m=joint.len();let n=joint[0].len();
    let px:Vec<f64>=(0..m).map(|i|joint[i].iter().sum()).collect();
    let py:Vec<f64>=(0..n).map(|j|joint.iter().map(|row|row[j]).sum()).collect();
    let mi=mutual_information(joint);let hx=entropy(&px);let hy=entropy(&py);
    if hx+hy==0.0{0.0}else{2.0*mi/(hx+hy)}
}

/// Perplexity = 2^H(X).
pub fn perplexity(probs:&[f64])->f64{exponential::exp2(entropy(probs))}

/// Bits needed to encode n equally likely outcomes.
pub fn hartley(n:u64)->f64{exponential::log2(n as f64)}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_entropy(){assert!(approx_eq(entropy(&[0.5,0.5]),1.0,1e-10));}
    #[test]fn test_entropy_certain(){assert!(approx_eq(entropy(&[1.0,0.0]),0.0,1e-10));}
    #[test]fn test_kl(){let p=[0.5,0.5];let q=[0.5,0.5];assert!(approx_eq(kl_divergence(&p,&q),0.0,1e-10));}
    #[test]fn test_cross_entropy(){let p=[0.5,0.5];assert!(approx_eq(cross_entropy(&p,&p),-0.5*exponential::ln(0.5)*2.0,1e-10));}
    #[test]fn test_hartley(){assert!(approx_eq(hartley(8),3.0,1e-10));}
}
