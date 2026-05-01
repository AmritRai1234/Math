//! Voting theory: electoral methods and social choice.

/// Plurality voting: most first-place votes wins. Returns winner index.
pub fn plurality(ballots:&[Vec<usize>],n_candidates:usize)->usize{
    let mut counts=vec![0usize;n_candidates];
    for b in ballots{if !b.is_empty(){counts[b[0]]+=1;}}
    counts.iter().enumerate().max_by_key(|(_,&c)|c).map(|(i,_)|i).unwrap_or(0)
}

/// Borda count: last place=0, first=n-1 points.
pub fn borda(ballots:&[Vec<usize>],n_candidates:usize)->Vec<f64>{
    let mut scores=vec![0.0;n_candidates];
    for b in ballots{let n=b.len();for(rank,&c) in b.iter().enumerate(){scores[c]+=(n-1-rank)as f64;}}scores
}

/// Borda winner.
pub fn borda_winner(ballots:&[Vec<usize>],n_candidates:usize)->usize{
    let scores=borda(ballots,n_candidates);
    scores.iter().enumerate().max_by(|(_,a),(_,b)|a.partial_cmp(b).unwrap()).map(|(i,_)|i).unwrap_or(0)
}

/// Condorcet winner (beats all others in pairwise). Returns Some(winner) or None.
pub fn condorcet(ballots:&[Vec<usize>],n_candidates:usize)->Option<usize>{
    let mut pairwise=vec![vec![0i64;n_candidates];n_candidates];
    for b in ballots{for i in 0..b.len(){for j in i+1..b.len(){pairwise[b[i]][b[j]]+=1;pairwise[b[j]][b[i]]-=1;}}}
    for c in 0..n_candidates{if(0..n_candidates).all(|o|o==c||pairwise[c][o]>0){return Some(c);}}None
}

/// Instant-runoff voting (ranked choice).
pub fn instant_runoff(ballots:&[Vec<usize>],n_candidates:usize)->usize{
    let mut eliminated=vec![false;n_candidates];let n_voters=ballots.len();
    loop{let mut counts=vec![0usize;n_candidates];
        for b in ballots{for &c in b{if !eliminated[c]{counts[c]+=1;break;}}}
        // Check majority
        for c in 0..n_candidates{if !eliminated[c]&&counts[c]*2>n_voters{return c;}}
        // Eliminate lowest
        let min_count=counts.iter().enumerate().filter(|(i,_)|!eliminated[*i]).map(|(_,&c)|c).min().unwrap_or(0);
        let mut found=false;for c in 0..n_candidates{if !eliminated[c]&&counts[c]==min_count{eliminated[c]=true;found=true;break;}}
        if !found{break;}
        // Check if only one remains
        let remaining:Vec<usize>=(0..n_candidates).filter(|&c|!eliminated[c]).collect();
        if remaining.len()<=1{return remaining.get(0).cloned().unwrap_or(0);}
    }
    0
}

/// Approval voting: each voter approves a subset.
pub fn approval(approvals:&[Vec<usize>],n_candidates:usize)->usize{
    let mut counts=vec![0usize;n_candidates];
    for ballot in approvals{for &c in ballot{counts[c]+=1;}}
    counts.iter().enumerate().max_by_key(|(_,&c)|c).map(|(i,_)|i).unwrap_or(0)
}

/// Copeland score: wins - losses in pairwise comparisons.
pub fn copeland(ballots:&[Vec<usize>],n_candidates:usize)->Vec<i64>{
    let mut pairwise=vec![vec![0i64;n_candidates];n_candidates];
    for b in ballots{for i in 0..b.len(){for j in i+1..b.len(){pairwise[b[i]][b[j]]+=1;}}}
    let mut scores=vec![0i64;n_candidates];
    for i in 0..n_candidates{for j in 0..n_candidates{if i!=j{
        if pairwise[i][j]>pairwise[j][i]{scores[i]+=1;}
        else if pairwise[i][j]<pairwise[j][i]{scores[i]-=1;}}}}
    scores
}

/// Score voting: each voter assigns scores. Returns totals.
pub fn score_voting(scores:&[Vec<f64>])->Vec<f64>{
    let n=scores[0].len();let mut totals=vec![0.0;n];
    for ballot in scores{for(i,&s) in ballot.iter().enumerate(){totals[i]+=s;}}totals
}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_plurality(){let b=vec![vec![0,1,2],vec![0,2,1],vec![1,0,2]];assert_eq!(plurality(&b,3),0);}
    #[test]fn test_borda(){let b=vec![vec![0,1,2],vec![0,1,2],vec![1,0,2]];let w=borda_winner(&b,3);assert_eq!(w,0);}
    #[test]fn test_condorcet(){let b=vec![vec![0,1,2],vec![0,1,2],vec![0,1,2]];assert_eq!(condorcet(&b,3),Some(0));}
    #[test]fn test_approval(){let a=vec![vec![0,1],vec![1,2],vec![0]];let w=approval(&a,3);assert!(w==0||w==1);}// 0 and 1 both have 2 approvals
}
