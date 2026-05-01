//! Data compression: Huffman coding, RLE, LZ77.

use std::collections::HashMap;

/// Run-length encoding.
pub fn rle_encode(data:&[u8])->Vec<(u8,u32)>{
    if data.is_empty(){return vec![];}
    let mut result=Vec::new();let mut cur=data[0];let mut count=1u32;
    for &b in &data[1..]{if b==cur{count+=1;}else{result.push((cur,count));cur=b;count=1;}}
    result.push((cur,count));result
}

/// Run-length decoding.
pub fn rle_decode(encoded:&[(u8,u32)])->Vec<u8>{
    encoded.iter().flat_map(|&(b,c)|(0..c).map(move|_|b)).collect()
}

/// Huffman tree node.
#[derive(Clone,Debug)]
enum HuffNode{Leaf(u8),Internal(Box<HuffNode>,Box<HuffNode>)}

/// Build Huffman codes from frequency table.
pub fn huffman_codes(data:&[u8])->HashMap<u8,String>{
    let mut freq=HashMap::new();for &b in data{*freq.entry(b).or_insert(0u32)+=1;}
    if freq.len()<=1{let mut m=HashMap::new();for(&k,_)in &freq{m.insert(k,"0".to_string());}return m;}
    // Build priority queue (simple sorted vec)
    let mut nodes:Vec<(u32,HuffNode)>=freq.iter().map(|(&b,&f)|(f,HuffNode::Leaf(b))).collect();
    while nodes.len()>1{nodes.sort_by(|a,b|b.0.cmp(&a.0));// sort descending so pop from end
        let(f1,n1)=nodes.pop().unwrap();let(f2,n2)=nodes.pop().unwrap();
        nodes.push((f1+f2,HuffNode::Internal(Box::new(n1),Box::new(n2))));}
    let mut codes=HashMap::new();
    fn traverse(node:&HuffNode,prefix:&str,codes:&mut HashMap<u8,String>){match node{
        HuffNode::Leaf(b)=>{codes.insert(*b,prefix.to_string());},
        HuffNode::Internal(l,r)=>{traverse(l,&format!("{}0",prefix),codes);traverse(r,&format!("{}1",prefix),codes);}}}
    if let Some((_,root))=nodes.pop(){traverse(&root,"",&mut codes);}codes
}

/// Compute compression ratio for Huffman encoding.
pub fn huffman_compression_ratio(data:&[u8])->f64{
    let codes=huffman_codes(data);let mut freq=HashMap::new();
    for &b in data{*freq.entry(b).or_insert(0u64)+=1;}
    let original_bits=data.len()*8;
    let compressed_bits:u64=freq.iter().map(|(b,&count)|count*codes[b].len()as u64).sum();
    if compressed_bits==0{return 1.0;}original_bits as f64/compressed_bits as f64
}

/// LZ77-style encoding: (offset, length, next_char).
pub fn lz77_encode(data:&[u8],window_size:usize)->Vec<(u16,u16,u8)>{
    let n=data.len();let mut i=0;let mut result=Vec::new();
    while i<n{let mut best_offset=0u16;let mut best_length=0u16;
        let start=if i>window_size{i-window_size}else{0};
        for j in start..i{let mut len=0;
            while i+len<n&&len<255&&data[j+len]==data[i+len]{len+=1;}
            if len>best_length as usize{best_length=len as u16;best_offset=(i-j)as u16;}}
        let next=if (i+best_length as usize)<n{data[i+best_length as usize]}else{0};
        result.push((best_offset,best_length,next));i+=best_length as usize+1;}result
}

/// LZ77 decode.
pub fn lz77_decode(encoded:&[(u16,u16,u8)])->Vec<u8>{
    let mut result=Vec::new();
    for &(offset,length,next) in encoded{let start=result.len()-offset as usize;
        for i in 0..length as usize{result.push(result[start+i]);}result.push(next);}result
}

/// Delta encoding.
pub fn delta_encode(data:&[i64])->Vec<i64>{
    if data.is_empty(){return vec![];}
    let mut result=vec![data[0]];for i in 1..data.len(){result.push(data[i]-data[i-1]);}result
}

/// Delta decoding.
pub fn delta_decode(encoded:&[i64])->Vec<i64>{
    if encoded.is_empty(){return vec![];}
    let mut result=vec![encoded[0]];for i in 1..encoded.len(){result.push(result[i-1]+encoded[i]);}result
}

/// Entropy (bits per symbol).
pub fn entropy(data:&[u8])->f64{
    let mut freq=HashMap::new();for &b in data{*freq.entry(b).or_insert(0u64)+=1;}
    let n=data.len()as f64;
    freq.values().map(|&c|{let p=c as f64/n;-p*crate::exponential::log2(p)}).sum()
}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_rle_roundtrip(){let data=b"aaabbbccdd";let enc=rle_encode(data);let dec=rle_decode(&enc);assert_eq!(dec,data);}
    #[test]fn test_huffman(){let codes=huffman_codes(b"aabbbcccc");assert!(codes[&b'c'].len()<=codes[&b'a'].len());}
    #[test]fn test_delta_roundtrip(){let data=vec![1,3,6,10,15];let enc=delta_encode(&data);let dec=delta_decode(&enc);assert_eq!(dec,data);}
    #[test]fn test_entropy(){let data:Vec<u8>=(0..256).map(|i|i as u8).collect();let e=entropy(&data);assert!((e-8.0).abs()<0.01);}
}
