//! Cryptographic math: modular arithmetic, RSA primitives, hashing.

use crate::number_theory::{mod_pow,mod_inv,is_prime};

/// Generate RSA key pair (small primes for demo). Returns ((e, n), (d, n)).
pub fn rsa_keygen(p:u64,q:u64)->((u64,u64),(u64,u64)){
    let n=p*q;let phi=(p-1)*(q-1);
    let e=65537u64;let d=mod_inv(e,phi).expect("e not invertible");
    ((e,n),(d,n))
}

/// RSA encrypt: c = m^e mod n.
pub fn rsa_encrypt(m:u64,e:u64,n:u64)->u64{mod_pow(m,e,n)}

/// RSA decrypt: m = c^d mod n.
pub fn rsa_decrypt(c:u64,d:u64,n:u64)->u64{mod_pow(c,d,n)}

/// Diffie-Hellman: compute public key g^a mod p.
pub fn dh_public(g:u64,private:u64,p:u64)->u64{mod_pow(g,private,p)}

/// Diffie-Hellman: compute shared secret.
pub fn dh_shared(other_public:u64,private:u64,p:u64)->u64{mod_pow(other_public,private,p)}

/// Baby-step giant-step discrete log: find x such that g^x ≡ h (mod p).
pub fn discrete_log(g:u64,h:u64,p:u64)->Option<u64>{
    let m=((p as f64).sqrt()as u64)+1;
    let mut table=std::collections::HashMap::new();
    let mut val=h;
    for j in 0..m{table.insert(val,j);val=(val as u128*g as u128%p as u128)as u64;}
    let gm=mod_pow(g,m,p);let mut gamma=1u64;
    for i in 1..=m{gamma=(gamma as u128*gm as u128%p as u128)as u64;
        if let Some(&j)=table.get(&gamma){let ans=i*m-j;if ans<p{return Some(ans);}}}
    None
}

/// FNV-1a hash (64-bit).
pub fn fnv1a_hash(data:&[u8])->u64{
    let mut hash=0xcbf29ce484222325u64;
    for &byte in data{hash^=byte as u64;hash=hash.wrapping_mul(0x100000001b3);}hash
}

/// Simple SHA-256-like mixing function (NOT cryptographically secure, demo only).
pub fn hash_mix(mut a:u32,mut b:u32,mut c:u32)->u32{
    a=a.wrapping_sub(c);a^=c.rotate_left(4);c=c.wrapping_add(b);
    b=b.wrapping_sub(a);b^=a.rotate_left(6);a=a.wrapping_add(c);
    c=c.wrapping_sub(b);c^=b.rotate_left(8);b=b.wrapping_add(a);
    a=a.wrapping_sub(c);a^=c.rotate_left(16);c=c.wrapping_add(b);
    c
}

/// Miller-Rabin witness test (exposed for crypto use).
pub fn is_probable_prime(n:u64,witnesses:&[u64])->bool{
    if n<2{return false;}if n==2||n==3{return true;}if n%2==0{return false;}
    let mut d=n-1;let mut r=0u32;while d%2==0{d/=2;r+=1;}
    'outer:for &a in witnesses{if a>=n{continue;}
        let mut x=mod_pow(a,d,n);if x==1||x==n-1{continue;}
        for _ in 0..r-1{x=mod_pow(x,2,n);if x==n-1{continue 'outer;}}
        return false;}
    true
}

/// Euler's criterion: check if a is a quadratic residue mod p.
pub fn is_quadratic_residue(a:u64,p:u64)->bool{mod_pow(a,(p-1)/2,p)==1}

/// Tonelli-Shanks: find sqrt(n) mod p. Returns None if n is not a QR.
pub fn mod_sqrt(n:u64,p:u64)->Option<u64>{
    if !is_quadratic_residue(n%p,p){return None;}
    if p%4==3{return Some(mod_pow(n,(p+1)/4,p));}
    // Tonelli-Shanks
    let mut q=p-1;let mut s=0u32;while q%2==0{q/=2;s+=1;}
    let mut z=2u64;while is_quadratic_residue(z,p){z+=1;}
    let mut m=s;let mut c=mod_pow(z,q,p);let mut t=mod_pow(n,q,p);let mut r=mod_pow(n,(q+1)/2,p);
    loop{if t==1{return Some(r);}
        let mut i=1u32;let mut tmp=mod_pow(t,2,p);while tmp!=1{tmp=mod_pow(tmp,2,p);i+=1;if i>=m{return None;}}
        let b=mod_pow(c,1<<(m-i-1),p);m=i;c=mod_pow(b,2,p);
        t=(t as u128*c as u128%p as u128)as u64;
        r=(r as u128*b as u128%p as u128)as u64;}
}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_rsa(){let((e,n),(d,_))=rsa_keygen(61,53);let msg=42u64;let c=rsa_encrypt(msg,e,n);let m=rsa_decrypt(c,d,n);assert_eq!(m,msg);}
    #[test]fn test_dh(){let(g,p)=(5u64,23u64);let a_priv=6;let b_priv=15;
        let a_pub=dh_public(g,a_priv,p);let b_pub=dh_public(g,b_priv,p);
        assert_eq!(dh_shared(b_pub,a_priv,p),dh_shared(a_pub,b_priv,p));}
    #[test]fn test_discrete_log(){let r=discrete_log(5,8,23);assert!(r.is_some());assert_eq!(mod_pow(5,r.unwrap(),23),8);}
    #[test]fn test_fnv(){let h1=fnv1a_hash(b"hello");let h2=fnv1a_hash(b"hello");assert_eq!(h1,h2);assert_ne!(fnv1a_hash(b"hello"),fnv1a_hash(b"world"));}
    #[test]fn test_mod_sqrt(){let r=mod_sqrt(2,7).unwrap();assert_eq!(r*r%7,2);}
}
