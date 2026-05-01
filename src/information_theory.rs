//! Information theory: entropy, KL divergence, mutual information, channel capacity.

use crate::exponential;use crate::arithmetic::abs;

/// Shannon entropy H(X) = -Σ p(x) log2 p(x).
pub fn entropy(probs:&[f64])->f64{
    probs.iter().filter(|&&p|p>0.0).map(|&p|-p*exponential::log2(p)).sum()
}

/// Joint entropy H(X,Y) from joint probability table.
pub fn joint_entropy(joint:&[Vec<f64>])->f64{
    joint.iter().flat_map(|row|row.iter()).filter(|&&p|p>0.0).map(|&p|-p*exponential::log2(p)).sum()
}

/// Conditional entropy H(X|Y) = H(X,Y) - H(Y).
pub fn conditional_entropy(joint:&[Vec<f64>])->f64{
    let hy=entropy(&marginal_y(joint));joint_entropy(joint)-hy
}

/// Mutual information I(X;Y) = H(X) + H(Y) - H(X,Y).
pub fn mutual_information(joint:&[Vec<f64>])->f64{
    let hx=entropy(&marginal_x(joint));let hy=entropy(&marginal_y(joint));
    hx+hy-joint_entropy(joint)
}

fn marginal_x(joint:&[Vec<f64>])->Vec<f64>{joint.iter().map(|row|row.iter().sum()).collect()}
fn marginal_y(joint:&[Vec<f64>])->Vec<f64>{
    let cols=joint[0].len();(0..cols).map(|j|joint.iter().map(|row|row[j]).sum()).collect()
}

/// KL divergence D_KL(P || Q) = Σ p(x) log(p(x)/q(x)).
pub fn kl_divergence(p:&[f64],q:&[f64])->f64{
    p.iter().zip(q.iter()).filter(|(&pi,&qi)|pi>0.0&&qi>0.0)
        .map(|(&pi,&qi)|pi*exponential::ln(pi/qi)).sum()
}

/// Symmetric KL divergence (Jensen-Shannon divergence).
pub fn js_divergence(p:&[f64],q:&[f64])->f64{
    let m:Vec<f64>=p.iter().zip(q.iter()).map(|(&pi,&qi)|(pi+qi)/2.0).collect();
    (kl_divergence(p,&m)+kl_divergence(q,&m))/2.0
}

/// Cross entropy H(P,Q) = -Σ p(x) log q(x).
pub fn cross_entropy(p:&[f64],q:&[f64])->f64{
    p.iter().zip(q.iter()).filter(|(&pi,&qi)|pi>0.0&&qi>0.0)
        .map(|(&pi,&qi)|-pi*exponential::log2(qi)).sum()
}

/// Binary entropy function H(p) = -p*log2(p) - (1-p)*log2(1-p).
pub fn binary_entropy(p:f64)->f64{
    if p<=0.0||p>=1.0{return 0.0;}-p*exponential::log2(p)-(1.0-p)*exponential::log2(1.0-p)
}

/// Channel capacity of binary symmetric channel: C = 1 - H(p).
pub fn bsc_capacity(p:f64)->f64{1.0-binary_entropy(p)}

/// Channel capacity of binary erasure channel: C = 1 - ε.
pub fn bec_capacity(epsilon:f64)->f64{1.0-epsilon}

/// Rényi entropy of order α.
pub fn renyi_entropy(probs:&[f64],alpha:f64)->f64{
    if(alpha-1.0).abs()<1e-10{return entropy(probs);}
    let sum:f64=probs.iter().map(|&p|if p>0.0{crate::powers::pow(p,alpha)}else{0.0}).sum();
    exponential::log2(sum)/(1.0-alpha)
}

/// Min-entropy (Rényi entropy α→∞).
pub fn min_entropy(probs:&[f64])->f64{
    -exponential::log2(probs.iter().cloned().fold(0.0f64,f64::max))
}

/// Hartley entropy (log of number of outcomes).
pub fn hartley_entropy(n:usize)->f64{exponential::log2(n as f64)}

/// Data processing inequality check: I(X;Z) ≤ I(X;Y) for Markov chain X→Y→Z.
pub fn data_processing_holds(ixy:f64,ixz:f64)->bool{ixz<=ixy+1e-10}

/// Normalized mutual information: NMI = 2*I(X;Y)/(H(X)+H(Y)).
pub fn normalized_mutual_info(joint:&[Vec<f64>])->f64{
    let mi=mutual_information(joint);let hx=entropy(&marginal_x(joint));let hy=entropy(&marginal_y(joint));
    if hx+hy<1e-30{return 0.0;}2.0*mi/(hx+hy)
}

/// Variation of information: VI = H(X|Y) + H(Y|X).
pub fn variation_of_information(joint:&[Vec<f64>])->f64{
    let hxy=joint_entropy(joint);let hx=entropy(&marginal_x(joint));let hy=entropy(&marginal_y(joint));
    2.0*hxy-hx-hy
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_fair_coin(){assert!(approx_eq(entropy(&[0.5,0.5]),1.0,1e-10));}
    #[test]fn test_certain(){assert!(approx_eq(entropy(&[1.0,0.0]),0.0,1e-10));}
    #[test]fn test_kl_same(){assert!(approx_eq(kl_divergence(&[0.5,0.5],&[0.5,0.5]),0.0,1e-10));}
    #[test]fn test_binary_entropy(){assert!(approx_eq(binary_entropy(0.5),1.0,1e-10));}
    #[test]fn test_bsc(){assert!(approx_eq(bsc_capacity(0.0),1.0,1e-10));}
    #[test]fn test_mutual_info(){let joint=vec![vec![0.25,0.25],vec![0.25,0.25]];// independent
        assert!(approx_eq(mutual_information(&joint),0.0,1e-10));}
}
