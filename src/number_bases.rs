//! Number base conversions, Roman numerals, balanced ternary.

/// Convert integer to string in given base (2-36).
pub fn to_base(mut n:u64,base:u64)->String{
    if n==0{return "0".to_string();}
    let digits=b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    let mut result=Vec::new();while n>0{result.push(digits[(n%base)as usize]);n/=base;}
    result.reverse();String::from_utf8(result).unwrap()
}

/// Convert string in given base to integer.
pub fn from_base(s:&str,base:u64)->u64{
    s.chars().fold(0u64,|acc,c|{let d=if c>='0'&&c<='9'{c as u64-'0' as u64}
        else if c>='A'&&c<='Z'{c as u64-'A' as u64+10}
        else if c>='a'&&c<='z'{c as u64-'a' as u64+10}else{0};acc*base+d})
}

/// Convert between arbitrary bases.
pub fn convert_base(s:&str,from:u64,to:u64)->String{to_base(from_base(s,from),to)}

/// Integer to binary string.
pub fn to_binary(n:u64)->String{to_base(n,2)}
/// Integer to octal string.
pub fn to_octal(n:u64)->String{to_base(n,8)}
/// Integer to hex string.
pub fn to_hex(n:u64)->String{to_base(n,16)}

/// Integer to Roman numerals.
pub fn to_roman(mut n:u32)->String{
    let vals=[(1000,"M"),(900,"CM"),(500,"D"),(400,"CD"),(100,"C"),(90,"XC"),
        (50,"L"),(40,"XL"),(10,"X"),(9,"IX"),(5,"V"),(4,"IV"),(1,"I")];
    let mut result=String::new();
    for &(val,sym) in &vals{while n>=val{result.push_str(sym);n-=val;}}result
}

/// Roman numerals to integer.
pub fn from_roman(s:&str)->u32{
    let val=|c|match c{'I'=>1,'V'=>5,'X'=>10,'L'=>50,'C'=>100,'D'=>500,'M'=>1000,_=>0};
    let chars:Vec<char>=s.chars().collect();let mut result=0u32;
    for i in 0..chars.len(){let v=val(chars[i]);
        if i+1<chars.len()&&v<val(chars[i+1]){result-=v;}else{result+=v;}}result
}

/// Integer to balanced ternary (-1, 0, 1 digits shown as T, 0, 1).
pub fn to_balanced_ternary(mut n:i64)->String{
    if n==0{return "0".to_string();}
    let neg=n<0;if neg{n=-n;}let mut digits=Vec::new();
    while n>0{let mut r=(n%3)as i8;n/=3;if r==2{r=-1;n+=1;}digits.push(r);}
    if neg{for d in &mut digits{*d=-*d;}}
    digits.reverse();digits.iter().map(|&d|match d{-1=>'T',0=>'0',1=>'1',_=>'?'}).collect()
}

/// Balanced ternary to integer.
pub fn from_balanced_ternary(s:&str)->i64{
    s.chars().fold(0i64,|acc,c|acc*3+match c{'T'=>-1,'0'=>0,'1'=>1,_=>0})
}

/// Gray code to binary.
pub fn gray_to_binary(gray:u64)->u64{crate::logic::from_gray(gray)}
/// Binary to Gray code.
pub fn binary_to_gray(n:u64)->u64{crate::logic::to_gray(n)}

/// BCD (Binary-Coded Decimal) encode.
pub fn to_bcd(n:u64)->u64{
    let mut result=0u64;let mut shift=0;let mut n=n;
    while n>0{result|=(n%10)<<shift;n/=10;shift+=4;}result
}

/// BCD decode.
pub fn from_bcd(bcd:u64)->u64{
    let mut result=0u64;let mut mult=1;let mut bcd=bcd;
    while bcd>0{result+=(bcd&0xF)*mult;bcd>>=4;mult*=10;}result
}

/// Two's complement representation.
pub fn twos_complement(n:i64,bits:u32)->u64{
    if n>=0{n as u64}else{let mask=(1u64<<bits)-1;((!(-n as u64))+1)&mask}
}

/// Factorial number system representation.
pub fn to_factoradic(mut n:u64)->Vec<u64>{
    let mut result=vec![0];let mut base=2;
    while n>0{result.push(n%base);n/=base;base+=1;}result.reverse();result
}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_base(){assert_eq!(to_base(255,16),"FF");assert_eq!(from_base("FF",16),255);}
    #[test]fn test_binary(){assert_eq!(to_binary(42),"101010");}
    #[test]fn test_roman(){assert_eq!(to_roman(1994),"MCMXCIV");assert_eq!(from_roman("MCMXCIV"),1994);}
    #[test]fn test_balanced_ternary(){let s=to_balanced_ternary(42);let back=from_balanced_ternary(&s);assert_eq!(back,42);}
    #[test]fn test_bcd(){assert_eq!(to_bcd(42),0x42);assert_eq!(from_bcd(0x42),42);}
    #[test]fn test_convert(){assert_eq!(convert_base("FF",16,2),"11111111");}
}
