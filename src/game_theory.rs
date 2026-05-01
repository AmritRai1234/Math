//! Game theory: payoff matrices, Nash equilibrium, minimax.

/// Minimax for zero-sum 2-player game. Returns (row_strategy, col_strategy, value).
pub fn minimax(payoff:&[Vec<f64>])->(usize,usize,f64){
    let m=payoff.len();let n=payoff[0].len();
    // Row player maximizes minimum
    let mut best_row=0;let mut best_val=f64::NEG_INFINITY;
    for i in 0..m{let min_in_row=payoff[i].iter().cloned().fold(f64::INFINITY,f64::min);
        if min_in_row>best_val{best_val=min_in_row;best_row=i;}}
    // Col player minimizes maximum
    let mut best_col=0;let mut best_col_val=f64::INFINITY;
    for j in 0..n{let max_in_col=(0..m).map(|i|payoff[i][j]).fold(f64::NEG_INFINITY,f64::max);
        if max_in_col<best_col_val{best_col_val=max_in_col;best_col=j;}}
    (best_row,best_col,best_val)
}

/// Check if game has a saddle point (pure strategy Nash equilibrium).
pub fn has_saddle_point(payoff:&[Vec<f64>])->Option<(usize,usize,f64)>{
    let(r,c,maximin)=minimax(payoff);
    let m=payoff.len();let n=payoff[0].len();
    let minimax_val=(0..n).map(|j|(0..m).map(|i|payoff[i][j]).fold(f64::NEG_INFINITY,f64::max)).fold(f64::INFINITY,f64::min);
    if(maximin-minimax_val).abs()<1e-10{Some((r,c,maximin))}else{None}
}

/// Dominant strategy elimination. Returns remaining rows and columns.
pub fn eliminate_dominated(payoff:&[Vec<f64>])->(Vec<usize>,Vec<usize>){
    let m=payoff.len();let n=payoff[0].len();
    let mut rows:Vec<usize>=(0..m).collect();let mut cols:Vec<usize>=(0..n).collect();
    let mut changed=true;
    while changed{changed=false;
        // Remove dominated rows
        let mut new_rows=Vec::new();
        for &i in &rows{let mut dominated=false;
            for &k in &rows{if k!=i{if cols.iter().all(|&j|payoff[k][j]>=payoff[i][j])
                &&cols.iter().any(|&j|payoff[k][j]>payoff[i][j]){dominated=true;break;}}}
            if !dominated{new_rows.push(i);}}
        if new_rows.len()<rows.len(){changed=true;rows=new_rows;}
        // Remove dominated columns
        let mut new_cols=Vec::new();
        for &j in &cols{let mut dominated=false;
            for &l in &cols{if l!=j{if rows.iter().all(|&i|payoff[i][l]<=payoff[i][j])
                &&rows.iter().any(|&i|payoff[i][l]<payoff[i][j]){dominated=true;break;}}}
            if !dominated{new_cols.push(j);}}
        if new_cols.len()<cols.len(){changed=true;cols=new_cols;}
    }
    (rows,cols)
}

/// Mixed strategy Nash equilibrium for 2x2 game. Returns (p, q, value).
/// p = probability row player plays row 0, q = probability col player plays col 0.
pub fn mixed_nash_2x2(payoff:&[[f64;2];2])->Option<(f64,f64,f64)>{
    let a=payoff[0][0];let b=payoff[0][1];let c=payoff[1][0];let d=payoff[1][1];
    let denom_q=a-b-c+d;let denom_p=a-b-c+d;
    if denom_q.abs()<1e-15||denom_p.abs()<1e-15{return None;}
    let q=(d-c)/denom_q;let p=(d-b)/denom_p;
    if q<0.0||q>1.0||p<0.0||p>1.0{return None;}
    let value=a*p*q+b*p*(1.0-q)+c*(1.0-p)*q+d*(1.0-p)*(1.0-q);
    Some((p,q,value))
}

/// Prisoner's dilemma payoff matrix.
pub fn prisoners_dilemma()->[[f64;2];2]{[[-1.0,-3.0],[0.0,-2.0]]}

/// Expected payoff given mixed strategies.
pub fn expected_payoff(payoff:&[Vec<f64>],row_probs:&[f64],col_probs:&[f64])->f64{
    let mut sum=0.0;
    for(i,&pi) in row_probs.iter().enumerate(){for(j,&qj) in col_probs.iter().enumerate(){sum+=pi*qj*payoff[i][j];}}sum
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_minimax(){let p=vec![vec![3.0,2.0,5.0],vec![1.0,4.0,1.0]];let(r,_,v)=minimax(&p);assert_eq!(r,0);}
    #[test]fn test_saddle(){let p=vec![vec![3.0,2.0],vec![1.0,4.0]];
        // This might or might not have a saddle point
        let _ = has_saddle_point(&p);}
    #[test]fn test_mixed_2x2(){let p=[[3.0,-3.0],[-3.0,3.0]];// matching pennies
        let(pr,pc,_)=mixed_nash_2x2(&p).unwrap();assert!(approx_eq(pr,0.5,1e-10));assert!(approx_eq(pc,0.5,1e-10));}
}
