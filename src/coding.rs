//! Error detection/correction codes: Hamming, CRC, parity.

/// Parity bit (even parity).
pub fn parity_bit(data:&[u8])->u8{let ones:u32=data.iter().map(|&b|b.count_ones()).sum();(ones%2)as u8}

/// Hamming(7,4) encode: 4 data bits → 7 bits.
pub fn hamming74_encode(data:u8)->u8{
    let d=[data&1,(data>>1)&1,(data>>2)&1,(data>>3)&1];
    let p1=d[0]^d[1]^d[3];let p2=d[0]^d[2]^d[3];let p3=d[1]^d[2]^d[3];
    p1|(p2<<1)|(d[0]<<2)|(p3<<3)|(d[1]<<4)|(d[2]<<5)|(d[3]<<6)
}

/// Hamming(7,4) decode with single-bit error correction. Returns (data, error_corrected).
pub fn hamming74_decode(mut code:u8)->(u8,bool){
    let bits:Vec<u8>=(0..7).map(|i|(code>>i)&1).collect();
    let s1=bits[0]^bits[2]^bits[4]^bits[6];let s2=bits[1]^bits[2]^bits[5]^bits[6];let s3=bits[3]^bits[4]^bits[5]^bits[6];
    let syndrome=s1|(s2<<1)|(s3<<2);
    let corrected=syndrome>0;if corrected&&(syndrome as usize)<=7{code^=1<<(syndrome-1);}
    let bits:Vec<u8>=(0..7).map(|i|(code>>i)&1).collect();
    let data=bits[2]|(bits[4]<<1)|(bits[5]<<2)|(bits[6]<<3);(data,corrected)
}

/// CRC-8 computation.
pub fn crc8(data:&[u8])->u8{
    let poly=0x07u8;let mut crc=0u8;
    for &byte in data{crc^=byte;for _ in 0..8{crc=if crc&0x80!=0{(crc<<1)^poly}else{crc<<1};}}crc
}

/// CRC-16 (CCITT).
pub fn crc16(data:&[u8])->u16{
    let poly=0x1021u16;let mut crc=0xFFFFu16;
    for &byte in data{crc^=(byte as u16)<<8;for _ in 0..8{crc=if crc&0x8000!=0{(crc<<1)^poly}else{crc<<1};}}crc
}

/// CRC-32 (Ethernet/ZIP).
pub fn crc32(data:&[u8])->u32{
    let poly=0xEDB88320u32;let mut crc=0xFFFFFFFFu32;
    for &byte in data{crc^=byte as u32;for _ in 0..8{crc=if crc&1!=0{(crc>>1)^poly}else{crc>>1};}}!crc
}

/// Luhn checksum (credit card validation).
pub fn luhn_check(digits:&[u8])->bool{
    let mut sum=0u32;let mut double=false;
    for &d in digits.iter().rev(){let mut n=d as u32;
        if double{n*=2;if n>9{n-=9;}}sum+=n;double=!double;}
    sum%10==0
}

/// Luhn compute check digit.
pub fn luhn_checkdigit(digits:&[u8])->u8{
    let mut extended=digits.to_vec();extended.push(0);
    for check in 0..10u8{extended[digits.len()]=check;if luhn_check(&extended){return check;}}0
}

/// ISBN-10 check digit.
pub fn isbn10_check(digits:&[u8])->bool{
    let sum:u32=digits.iter().enumerate().map(|(i,&d)|d as u32*(10-i as u32)).sum();sum%11==0
}

/// Repetition code encode (triple repetition).
pub fn rep3_encode(bit:u8)->u8{bit&1|((bit&1)<<1)|((bit&1)<<2)}

/// Repetition code decode (majority vote).
pub fn rep3_decode(code:u8)->u8{let ones=(code&1)+((code>>1)&1)+((code>>2)&1);if ones>=2{1}else{0}}

/// Hamming distance between two byte slices.
pub fn byte_hamming_distance(a:&[u8],b:&[u8])->u32{
    a.iter().zip(b.iter()).map(|(&x,&y)|(x^y).count_ones()).sum()
}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_hamming74(){let encoded=hamming74_encode(0b1011);let(decoded,_)=hamming74_decode(encoded);assert_eq!(decoded,0b1011);}
    #[test]fn test_hamming74_correct(){let encoded=hamming74_encode(0b1011);let corrupted=encoded^0b0100;// flip bit 2
        let(decoded,corrected)=hamming74_decode(corrupted);assert_eq!(decoded,0b1011);assert!(corrected);}
    #[test]fn test_crc32(){let crc=crc32(b"123456789");assert_eq!(crc,0xCBF43926);}
    #[test]fn test_luhn(){assert!(luhn_check(&[4,5,3,9,1,4,8,8,0,3,4,3,6,4,6,7]));}
    #[test]fn test_rep3(){assert_eq!(rep3_decode(rep3_encode(1)),1);assert_eq!(rep3_decode(0b110),1);}
    #[test]fn test_parity(){assert_eq!(parity_bit(&[0b10110101]),1);}
}
