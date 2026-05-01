//! Boolean algebra and propositional logic.

/// Truth table for a function of n boolean variables.
pub fn truth_table(n:u32,f:&dyn Fn(&[bool])->bool)->Vec<(Vec<bool>,bool)>{
    let rows=1u32<<n;let mut table=Vec::with_capacity(rows as usize);
    for i in 0..rows{let inputs:Vec<bool>=(0..n).map(|j|(i>>(n-1-j))&1==1).collect();
        let output=f(&inputs);table.push((inputs,output));}table
}

/// AND.
pub fn and(a:bool,b:bool)->bool{a&&b}
/// OR.
pub fn or(a:bool,b:bool)->bool{a||b}
/// NOT.
pub fn not(a:bool)->bool{!a}
/// XOR.
pub fn xor(a:bool,b:bool)->bool{a^b}
/// NAND.
pub fn nand(a:bool,b:bool)->bool{!(a&&b)}
/// NOR.
pub fn nor(a:bool,b:bool)->bool{!(a||b)}
/// XNOR.
pub fn xnor(a:bool,b:bool)->bool{!(a^b)}
/// Implication: a → b = ¬a ∨ b.
pub fn implies(a:bool,b:bool)->bool{!a||b}
/// Biconditional: a ↔ b.
pub fn iff(a:bool,b:bool)->bool{a==b}

/// Majority function: true if more than half are true.
pub fn majority(inputs:&[bool])->bool{
    let count=inputs.iter().filter(|&&x|x).count();count*2>inputs.len()
}

/// Parity: true if odd number of trues.
pub fn parity(inputs:&[bool])->bool{
    inputs.iter().filter(|&&x|x).count()%2==1
}

/// Convert boolean function to minterm (SOP) canonical form.
/// Returns vec of minterms (input combinations where output is true).
pub fn minterms(n:u32,f:&dyn Fn(&[bool])->bool)->Vec<u32>{
    let table=truth_table(n,f);
    table.iter().enumerate().filter(|(_,(_,out))|*out).map(|(i,_)|i as u32).collect()
}

/// Convert boolean function to maxterm (POS) canonical form.
pub fn maxterms(n:u32,f:&dyn Fn(&[bool])->bool)->Vec<u32>{
    let table=truth_table(n,f);
    table.iter().enumerate().filter(|(_,(_,out))|!*out).map(|(i,_)|i as u32).collect()
}

/// Check if function is tautology (always true).
pub fn is_tautology(n:u32,f:&dyn Fn(&[bool])->bool)->bool{
    truth_table(n,f).iter().all(|(_,out)|*out)
}

/// Check if function is contradiction (always false).
pub fn is_contradiction(n:u32,f:&dyn Fn(&[bool])->bool)->bool{
    truth_table(n,f).iter().all(|(_,out)|!*out)
}

/// Check if function is satisfiable.
pub fn is_satisfiable(n:u32,f:&dyn Fn(&[bool])->bool)->bool{
    truth_table(n,f).iter().any(|(_,out)|*out)
}

/// Check if two functions are logically equivalent.
pub fn equivalent(n:u32,f:&dyn Fn(&[bool])->bool,g:&dyn Fn(&[bool])->bool)->bool{
    let tf=truth_table(n,f);let tg=truth_table(n,g);
    tf.iter().zip(tg.iter()).all(|((_,a),(_,b))|a==b)
}

/// Hamming weight (number of 1s in binary).
pub fn hamming_weight(n:u64)->u32{crate::bits::popcount(n)}

/// Hamming distance between two values.
pub fn hamming_distance(a:u64,b:u64)->u32{crate::bits::popcount(a^b)}

/// Gray code.
pub fn to_gray(n:u64)->u64{n^(n>>1)}
pub fn from_gray(mut g:u64)->u64{let mut n=g;let mut mask=g>>1;while mask!=0{n^=mask;mask>>=1;}n}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_gates(){assert!(and(true,true));assert!(!and(true,false));assert!(or(true,false));assert!(xor(true,false));assert!(!xor(true,true));}
    #[test]fn test_truth_table(){let t=truth_table(2,&|v|and(v[0],v[1]));assert_eq!(t.len(),4);assert_eq!(t[3].1,true);}
    #[test]fn test_tautology(){assert!(is_tautology(1,&|v|or(v[0],not(v[0]))));}
    #[test]fn test_gray(){for i in 0..16{assert_eq!(from_gray(to_gray(i)),i);}}
    #[test]fn test_hamming(){assert_eq!(hamming_distance(0b1011,0b1101),2);}
}
