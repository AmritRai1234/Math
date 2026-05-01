//! String metrics: edit distances, similarity measures.

/// Levenshtein distance (edit distance).
pub fn levenshtein(a:&str,b:&str)->usize{
    let(m,n)=(a.len(),b.len());let a:Vec<char>=a.chars().collect();let b:Vec<char>=b.chars().collect();
    let mut dp=vec![vec![0usize;n+1];m+1];
    for i in 0..=m{dp[i][0]=i;}for j in 0..=n{dp[0][j]=j;}
    for i in 1..=m{for j in 1..=n{let cost=if a[i-1]==b[j-1]{0}else{1};
        dp[i][j]=(dp[i-1][j]+1).min(dp[i][j-1]+1).min(dp[i-1][j-1]+cost);}}dp[m][n]
}

/// Damerau-Levenshtein distance (includes transpositions).
pub fn damerau_levenshtein(a:&str,b:&str)->usize{
    let(m,n)=(a.len(),b.len());let a:Vec<char>=a.chars().collect();let b:Vec<char>=b.chars().collect();
    let mut dp=vec![vec![0usize;n+1];m+1];
    for i in 0..=m{dp[i][0]=i;}for j in 0..=n{dp[0][j]=j;}
    for i in 1..=m{for j in 1..=n{let cost=if a[i-1]==b[j-1]{0}else{1};
        dp[i][j]=(dp[i-1][j]+1).min(dp[i][j-1]+1).min(dp[i-1][j-1]+cost);
        if i>1&&j>1&&a[i-1]==b[j-2]&&a[i-2]==b[j-1]{dp[i][j]=dp[i][j].min(dp[i-2][j-2]+cost);}}}dp[m][n]
}

/// Hamming distance (for equal-length strings).
pub fn hamming(a:&str,b:&str)->usize{a.chars().zip(b.chars()).filter(|&(x,y)|x!=y).count()}

/// Longest common subsequence length.
pub fn lcs_length(a:&str,b:&str)->usize{
    let(m,n)=(a.len(),b.len());let a:Vec<char>=a.chars().collect();let b:Vec<char>=b.chars().collect();
    let mut dp=vec![vec![0usize;n+1];m+1];
    for i in 1..=m{for j in 1..=n{dp[i][j]=if a[i-1]==b[j-1]{dp[i-1][j-1]+1}else{dp[i-1][j].max(dp[i][j-1])};}}dp[m][n]
}

/// Longest common subsequence string.
pub fn lcs(a:&str,b:&str)->String{
    let(m,n)=(a.len(),b.len());let a:Vec<char>=a.chars().collect();let b:Vec<char>=b.chars().collect();
    let mut dp=vec![vec![0usize;n+1];m+1];
    for i in 1..=m{for j in 1..=n{dp[i][j]=if a[i-1]==b[j-1]{dp[i-1][j-1]+1}else{dp[i-1][j].max(dp[i][j-1])};}}
    let mut result=Vec::new();let(mut i,mut j)=(m,n);
    while i>0&&j>0{if a[i-1]==b[j-1]{result.push(a[i-1]);i-=1;j-=1;}else if dp[i-1][j]>dp[i][j-1]{i-=1;}else{j-=1;}}
    result.reverse();result.iter().collect()
}

/// Jaccard similarity (on character sets).
pub fn jaccard(a:&str,b:&str)->f64{
    let sa:std::collections::HashSet<char>=a.chars().collect();
    let sb:std::collections::HashSet<char>=b.chars().collect();
    let inter=sa.intersection(&sb).count()as f64;let union=sa.union(&sb).count()as f64;
    if union==0.0{1.0}else{inter/union}
}

/// Cosine similarity (on character frequency vectors).
pub fn cosine_similarity(a:&str,b:&str)->f64{
    let mut freq_a=std::collections::HashMap::new();let mut freq_b=std::collections::HashMap::new();
    for c in a.chars(){*freq_a.entry(c).or_insert(0i64)+=1;}
    for c in b.chars(){*freq_b.entry(c).or_insert(0i64)+=1;}
    let mut dot=0i64;let mut mag_a=0i64;let mut mag_b=0i64;
    let all_chars:std::collections::HashSet<char>=freq_a.keys().chain(freq_b.keys()).cloned().collect();
    for c in all_chars{let fa=*freq_a.get(&c).unwrap_or(&0);let fb=*freq_b.get(&c).unwrap_or(&0);
        dot+=fa*fb;mag_a+=fa*fa;mag_b+=fb*fb;}
    if mag_a==0||mag_b==0{return 0.0;}
    dot as f64/crate::powers::sqrt(mag_a as f64*mag_b as f64)
}

/// Jaro similarity.
pub fn jaro(a:&str,b:&str)->f64{
    let(a,b)=(a.chars().collect::<Vec<_>>(),b.chars().collect::<Vec<_>>());
    let(la,lb)=(a.len(),b.len());if la==0&&lb==0{return 1.0;}if la==0||lb==0{return 0.0;}
    let match_dist=(la.max(lb)/2).saturating_sub(1);
    let mut a_matches=vec![false;la];let mut b_matches=vec![false;lb];let mut matches=0.0;
    for i in 0..la{let start=if i>match_dist{i-match_dist}else{0};let end=(i+match_dist+1).min(lb);
        for j in start..end{if !b_matches[j]&&a[i]==b[j]{a_matches[i]=true;b_matches[j]=true;matches+=1.0;break;}}}
    if matches==0.0{return 0.0;}
    let mut transpositions=0.0;let mut k=0;
    for i in 0..la{if a_matches[i]{while!b_matches[k]{k+=1;}if a[i]!=b[k]{transpositions+=1.0;}k+=1;}}
    (matches/la as f64+matches/lb as f64+(matches-transpositions/2.0)/matches)/3.0
}

/// Jaro-Winkler similarity.
pub fn jaro_winkler(a:&str,b:&str)->f64{
    let j=jaro(a,b);let prefix=a.chars().zip(b.chars()).take(4).take_while(|&(x,y)|x==y).count();
    j+prefix as f64*0.1*(1.0-j)
}

/// Normalized Levenshtein similarity (0 to 1).
pub fn levenshtein_similarity(a:&str,b:&str)->f64{
    let max_len=a.len().max(b.len());if max_len==0{return 1.0;}1.0-levenshtein(a,b)as f64/max_len as f64
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_levenshtein(){assert_eq!(levenshtein("kitten","sitting"),3);}
    #[test]fn test_hamming(){assert_eq!(hamming("karolin","kathrin"),3);}
    #[test]fn test_lcs(){assert_eq!(lcs_length("ABCBDAB","BDCAB"),4);} // multiple valid LCS exist
    #[test]fn test_jaccard(){assert!(approx_eq(jaccard("abc","bcd"),0.5,1e-10));}
    #[test]fn test_jaro(){let j=jaro("martha","marhta");assert!(j>0.94);}
    #[test]fn test_identity(){assert!(approx_eq(levenshtein_similarity("hello","hello"),1.0,1e-10));}
}
