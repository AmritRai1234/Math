//! Finite fields: GF(p) arithmetic, GF(2^n) polynomial arithmetic.
fn mod_inv_internal(a:u64,p:u64)->u64{
    let(mut old_r,mut r)=(a as i64,p as i64);let(mut old_s,mut s)=(1i64,0i64);
    while r!=0{let q=old_r/r;let tmp=r;r=old_r-q*r;old_r=tmp;let tmp=s;s=old_s-q*s;old_s=tmp;}
    ((old_s%p as i64+p as i64)%p as i64)as u64
}


/// GF(p) element: arithmetic modulo prime p.
#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub struct GFp{pub val:u64,pub p:u64}

impl GFp{
    pub fn new(val:u64,p:u64)->Self{Self{val:val%p,p}}
    pub fn add(self,other:Self)->Self{Self{val:(self.val+other.val)%self.p,p:self.p}}
    pub fn sub(self,other:Self)->Self{Self{val:(self.val+self.p-other.val)%self.p,p:self.p}}
    pub fn mul(self,other:Self)->Self{Self{val:(self.val*other.val)%self.p,p:self.p}}
    pub fn inv(self)->Self{Self{val:mod_inv_internal(self.val,self.p),p:self.p}}
    pub fn div(self,other:Self)->Self{self.mul(other.inv())}
    pub fn pow(self,exp:u64)->Self{Self{val:crate::number_theory::mod_pow(self.val,exp,self.p),p:self.p}}
    pub fn neg(self)->Self{Self{val:(self.p-self.val)%self.p,p:self.p}}
    pub fn is_zero(self)->bool{self.val==0}
}

/// GF(2) polynomial: represented as bits of u64 (coefficient of x^i is bit i).
#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub struct GF2Poly(pub u64);

impl GF2Poly{
    pub fn new(bits:u64)->Self{Self(bits)}
    pub fn degree(self)->u32{if self.0==0{0}else{63-self.0.leading_zeros()}}

    /// Addition in GF(2)[x]: XOR.
    pub fn add(self,other:Self)->Self{Self(self.0^other.0)}
    /// Subtraction = addition in GF(2).
    pub fn sub(self,other:Self)->Self{self.add(other)}

    /// Multiplication in GF(2)[x].
    pub fn mul(self,other:Self)->Self{
        let mut result=0u64;let mut b=other.0;let mut shift=0;
        while b>0{if b&1==1{result^=self.0<<shift;}b>>=1;shift+=1;}Self(result)
    }

    /// Division: returns (quotient, remainder).
    pub fn divmod(self,other:Self)->(Self,Self){
        if other.0==0{panic!("division by zero polynomial");}
        let mut remainder=self.0;let divisor=other.0;let dd=other.degree();
        let mut quotient=0u64;
        while remainder!=0&&(63-remainder.leading_zeros())>=dd{
            let shift=(63-remainder.leading_zeros())-dd;
            quotient^=1<<shift;remainder^=divisor<<shift;}
        (Self(quotient),Self(remainder))
    }

    /// Modular reduction.
    pub fn modulo(self,modulus:Self)->Self{self.divmod(modulus).1}

    /// Multiplication modulo irreducible polynomial.
    pub fn mul_mod(self,other:Self,modulus:Self)->Self{self.mul(other).modulo(modulus)}

    /// GCD of two GF(2) polynomials.
    pub fn gcd(mut a:Self,mut b:Self)->Self{while b.0!=0{let r=a.divmod(b).1;a=b;b=r;}a}

    /// Check if polynomial is irreducible (brute force for small degree).
    pub fn is_irreducible(self)->bool{
        let d=self.degree();if d<=1{return d==1;}
        for i in 2..=(1u64<<(d/2)){let trial=Self(i);
            if trial.degree()>=1&&trial.degree()<d{if self.divmod(trial).1.0==0{return false;}}}true
    }
}

/// GF(2^n) element: polynomial modulo irreducible.
#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub struct GF2n{pub poly:GF2Poly,pub modulus:GF2Poly}

impl GF2n{
    pub fn new(val:u64,modulus:u64)->Self{
        let m=GF2Poly(modulus);Self{poly:GF2Poly(val).modulo(m),modulus:m}
    }
    pub fn add(self,other:Self)->Self{Self{poly:self.poly.add(other.poly),modulus:self.modulus}}
    pub fn mul(self,other:Self)->Self{Self{poly:self.poly.mul_mod(other.poly,self.modulus),modulus:self.modulus}}
    pub fn pow(self,mut exp:u64)->Self{
        let mut result=Self{poly:GF2Poly(1),modulus:self.modulus};let mut base=self;
        while exp>0{if exp&1==1{result=result.mul(base);}base=base.mul(base);exp>>=1;}result
    }
}

/// Find an irreducible polynomial of degree n over GF(2).
pub fn find_irreducible(n:u32)->GF2Poly{
    let start=1u64<<n;let end=start<<1;
    for i in start..end{let p=GF2Poly(i|1);// must have constant term
        if p.is_irreducible(){return p;}}
    GF2Poly(start|1) // fallback
}

/// Generate all elements of GF(p).
pub fn gf_elements(p:u64)->Vec<GFp>{(0..p).map(|v|GFp::new(v,p)).collect()}

/// Multiplicative order of element in GF(p).
pub fn multiplicative_order(a:GFp)->u64{
    if a.is_zero(){return 0;}let mut x=a;
    for i in 1..a.p{if x.val==1{return i;}x=x.mul(a);}a.p-1
}

/// Check if element is primitive root (generator) of GF(p)*.
pub fn is_primitive_root(a:GFp)->bool{multiplicative_order(a)==a.p-1}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_gfp(){let a=GFp::new(3,7);let b=GFp::new(5,7);assert_eq!(a.add(b).val,1);assert_eq!(a.mul(b).val,1);}
    #[test]fn test_gfp_inv(){let a=GFp::new(3,7);assert_eq!(a.mul(a.inv()).val,1);}
    #[test]fn test_gf2poly_mul(){let a=GF2Poly(0b111);let b=GF2Poly(0b11);let c=a.mul(b);assert_eq!(c.0,0b1001);}
    #[test]fn test_irreducible(){assert!(GF2Poly(0b111).is_irreducible());} // x²+x+1
    #[test]fn test_primitive(){assert!(is_primitive_root(GFp::new(3,7)));} // 3 is primitive root mod 7
}
