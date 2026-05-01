//! Bit manipulation functions.

/// Count set bits (population count).
pub fn popcount(mut n:u64)->u32{let mut c=0;while n!=0{c+=1;n&=n-1;}c}

/// Count leading zeros.
pub fn clz(n:u64)->u32{if n==0{return 64;}let mut c=0u32;let mut x=n;while x&(1u64<<63)==0{c+=1;x<<=1;}c}

/// Count trailing zeros.
pub fn ctz(n:u64)->u32{if n==0{return 64;}let mut c=0u32;let mut x=n;while x&1==0{c+=1;x>>=1;}c}

/// Next power of 2 >= n.
pub fn next_pow2(n:u64)->u64{if n==0{return 1;}let mut v=n-1;v|=v>>1;v|=v>>2;v|=v>>4;v|=v>>8;v|=v>>16;v|=v>>32;v+1}

/// Check if n is a power of 2.
pub fn is_pow2(n:u64)->bool{n!=0&&(n&(n-1))==0}

/// Rotate bits left.
pub fn rotate_left(n:u64,k:u32)->u64{let k=k%64;(n<<k)|(n>>(64-k))}

/// Rotate bits right.
pub fn rotate_right(n:u64,k:u32)->u64{let k=k%64;(n>>k)|(n<<(64-k))}

/// Reverse bit order.
pub fn reverse_bits(mut n:u64)->u64{let mut r=0u64;for _ in 0..64{r=(r<<1)|(n&1);n>>=1;}r}

/// Bit parity (1 if odd number of set bits, 0 if even).
pub fn parity(n:u64)->u32{popcount(n)%2}

/// Extract bits [low, high] inclusive.
pub fn extract_bits(n:u64,low:u32,high:u32)->u64{let mask=((1u64<<(high-low+1))-1)<<low;(n&mask)>>low}

/// Set bit at position.
pub fn set_bit(n:u64,pos:u32)->u64{n|(1u64<<pos)}

/// Clear bit at position.
pub fn clear_bit(n:u64,pos:u32)->u64{n&!(1u64<<pos)}

/// Toggle bit at position.
pub fn toggle_bit(n:u64,pos:u32)->u64{n^(1u64<<pos)}

/// Check if bit is set.
pub fn test_bit(n:u64,pos:u32)->bool{(n>>pos)&1==1}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_popcount(){assert_eq!(popcount(0b10110),3);assert_eq!(popcount(0),0);}
    #[test]fn test_clz(){assert_eq!(clz(1),63);assert_eq!(clz(0),64);}
    #[test]fn test_ctz(){assert_eq!(ctz(8),3);assert_eq!(ctz(0),64);}
    #[test]fn test_next_pow2(){assert_eq!(next_pow2(5),8);assert_eq!(next_pow2(8),8);assert_eq!(next_pow2(1),1);}
    #[test]fn test_is_pow2(){assert!(is_pow2(8));assert!(!is_pow2(6));}
    #[test]fn test_reverse(){assert_eq!(reverse_bits(1),1u64<<63);}
    #[test]fn test_parity(){assert_eq!(parity(0b111),1);assert_eq!(parity(0b110),0);}
}
