//! Big integer arithmetic: arbitrary-precision add, multiply, factorial.

/// Big integer stored as little-endian digits (base 10^9).
#[derive(Clone,Debug,PartialEq,Eq)]
pub struct BigInt{pub digits:Vec<u64>,pub negative:bool}
const BASE:u64=1_000_000_000;

impl BigInt{
    pub fn zero()->Self{Self{digits:vec![0],negative:false}}
    pub fn one()->Self{Self{digits:vec![1],negative:false}}
    pub fn from_u64(n:u64)->Self{if n==0{Self::zero()}else if n<BASE{Self{digits:vec![n],negative:false}}
        else{Self{digits:vec![n%BASE,n/BASE],negative:false}}}
    pub fn from_str(s:&str)->Self{
        let(negative,s)=if s.starts_with('-'){(true,&s[1..])}else{(false,s)};
        let mut digits=Vec::new();let bytes=s.as_bytes();let n=bytes.len();
        let mut i=n;while i>0{let start=if i>=9{i-9}else{0};
            let chunk=&s[start..i];digits.push(chunk.parse::<u64>().unwrap_or(0));i=start;}
        while digits.len()>1&&*digits.last().unwrap()==0{digits.pop();}
        Self{digits,negative}
    }
    pub fn to_string(&self)->String{
        if self.digits.is_empty()||self.digits==[0]{return "0".into();}
        let mut s=if self.negative{"-".to_string()}else{String::new()};
        s+=&self.digits.last().unwrap().to_string();
        for d in self.digits.iter().rev().skip(1){s+=&format!("{:09}",d);}s
    }
    fn trim(&mut self){while self.digits.len()>1&&*self.digits.last().unwrap()==0{self.digits.pop();}}
}

/// Add two positive BigInts.
pub fn bigint_add(a:&BigInt,b:&BigInt)->BigInt{
    let n=a.digits.len().max(b.digits.len());let mut result=Vec::with_capacity(n+1);let mut carry=0u64;
    for i in 0..n{let da=if i<a.digits.len(){a.digits[i]}else{0};
        let db=if i<b.digits.len(){b.digits[i]}else{0};
        let sum=da+db+carry;result.push(sum%BASE);carry=sum/BASE;}
    if carry>0{result.push(carry);}let mut r=BigInt{digits:result,negative:false};r.trim();r
}

/// Multiply two BigInts.
pub fn bigint_mul(a:&BigInt,b:&BigInt)->BigInt{
    let n=a.digits.len()+b.digits.len();let mut result=vec![0u64;n];
    for i in 0..a.digits.len(){let mut carry=0u64;for j in 0..b.digits.len(){
        let prod=a.digits[i] as u128*b.digits[j]as u128+result[i+j]as u128+carry as u128;
        result[i+j]=(prod%BASE as u128)as u64;carry=(prod/BASE as u128)as u64;}
        if carry>0{result[i+b.digits.len()]+=carry;}}
    let mut r=BigInt{digits:result,negative:a.negative!=b.negative};r.trim();r
}

/// BigInt factorial.
pub fn bigint_factorial(n:u64)->BigInt{
    let mut result=BigInt::one();for i in 2..=n{result=bigint_mul(&result,&BigInt::from_u64(i));}result
}

/// BigInt power.
pub fn bigint_pow(base:&BigInt,mut exp:u64)->BigInt{
    let mut result=BigInt::one();let mut b=base.clone();
    while exp>0{if exp%2==1{result=bigint_mul(&result,&b);}b=bigint_mul(&b,&b);exp/=2;}result
}

/// Compare BigInts (positive only).
pub fn bigint_cmp(a:&BigInt,b:&BigInt)->std::cmp::Ordering{
    use std::cmp::Ordering::*;
    if a.digits.len()!=b.digits.len(){return a.digits.len().cmp(&b.digits.len());}
    for i in (0..a.digits.len()).rev(){if a.digits[i]!=b.digits[i]{return a.digits[i].cmp(&b.digits[i]);}}Equal
}

/// Subtract (a - b), assuming a >= b.
pub fn bigint_sub(a:&BigInt,b:&BigInt)->BigInt{
    let mut result=Vec::with_capacity(a.digits.len());let mut borrow=0i64;
    for i in 0..a.digits.len(){let da=a.digits[i]as i64;
        let db=if i<b.digits.len(){b.digits[i]as i64}else{0};
        let diff=da-db-borrow;if diff<0{result.push((diff+BASE as i64)as u64);borrow=1;}
        else{result.push(diff as u64);borrow=0;}}
    let mut r=BigInt{digits:result,negative:false};r.trim();r
}

/// Number of digits.
pub fn bigint_num_digits(a:&BigInt)->usize{
    if a.digits.is_empty(){return 0;}
    let last=a.digits.last().unwrap().to_string().len();
    (a.digits.len()-1)*9+last
}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_add(){let a=BigInt::from_u64(999999999);let b=BigInt::from_u64(1);
        let c=bigint_add(&a,&b);assert_eq!(c.to_string(),"1000000000");}
    #[test]fn test_mul(){let a=BigInt::from_u64(12345);let b=BigInt::from_u64(67890);
        let c=bigint_mul(&a,&b);assert_eq!(c.to_string(),"838102050");}
    #[test]fn test_factorial(){let f=bigint_factorial(20);assert_eq!(f.to_string(),"2432902008176640000");}
    #[test]fn test_pow(){let b=BigInt::from_u64(2);let r=bigint_pow(&b,10);assert_eq!(r.to_string(),"1024");}
}
