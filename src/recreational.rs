//! Recreational math: magic squares, Sudoku solver, knight's tour, Collatz.

/// Generate odd-order magic square (Siamese method).
pub fn magic_square(n:usize)->Vec<Vec<u32>>{
    assert!(n%2==1&&n>=3);let mut sq=vec![vec![0u32;n];n];
    let(mut i,mut j)=(0,n/2);
    for num in 1..=n*n{sq[i][j]=num as u32;
        let(ni,nj)=((i+n-1)%n,(j+1)%n);
        if sq[ni][nj]!=0{i=(i+1)%n;}else{i=ni;j=nj;}}sq
}

/// Magic constant for n×n magic square.
pub fn magic_constant(n:u32)->u32{n*(n*n+1)/2}

/// Verify magic square.
pub fn is_magic_square(sq:&[Vec<u32>])->bool{
    let n=sq.len();let mc=magic_constant(n as u32);
    // Check rows
    for row in sq{if row.iter().sum::<u32>()!=mc{return false;}}
    // Check columns
    for j in 0..n{if (0..n).map(|i|sq[i][j]).sum::<u32>()!=mc{return false;}}
    // Check diagonals
    if (0..n).map(|i|sq[i][i]).sum::<u32>()!=mc{return false;}
    if (0..n).map(|i|sq[i][n-1-i]).sum::<u32>()!=mc{return false;}
    true
}

/// Collatz sequence.
pub fn collatz(mut n:u64)->Vec<u64>{
    let mut seq=vec![n];while n!=1{n=if n%2==0{n/2}else{3*n+1};seq.push(n);}seq
}

/// Collatz stopping time.
pub fn collatz_stopping_time(n:u64)->u64{collatz(n).len()as u64-1}

/// Happy number check.
pub fn is_happy(mut n:u64)->bool{
    let mut seen=std::collections::HashSet::new();
    while n!=1{if !seen.insert(n){return false;}
        let mut sum=0;while n>0{let d=n%10;sum+=d*d;n/=10;}n=sum;}true
}

/// Perfect number check: n equals sum of its proper divisors.
pub fn is_perfect(n:u64)->bool{if n<2{return false;}
    let sum:u64=(1..n).filter(|&i|n%i==0).sum();sum==n
}

/// Amicable pair check.
pub fn are_amicable(a:u64,b:u64)->bool{
    let sum_div=|n:u64|(1..n).filter(|&i|n%i==0).sum::<u64>();
    a!=b&&sum_div(a)==b&&sum_div(b)==a
}

/// Narcissistic number check (Armstrong number).
pub fn is_narcissistic(n:u64)->bool{
    let s=format!("{}",n);let d=s.len()as u32;
    let sum:u64=s.chars().map(|c|(c as u64-'0' as u64).pow(d)).sum();sum==n
}

/// Kaprekar's routine: sort digits desc - asc, repeat.
pub fn kaprekar_routine(mut n:u32,steps:usize)->Vec<u32>{
    let mut seq=vec![n];
    for _ in 0..steps{let s=format!("{:04}",n);let mut digits:Vec<u8>=s.bytes().collect();
        digits.sort();let asc:u32=digits.iter().fold(0u32,|acc,&d|acc*10+(d-b'0')as u32);
        digits.reverse();let desc:u32=digits.iter().fold(0u32,|acc,&d|acc*10+(d-b'0')as u32);
        n=desc-asc;seq.push(n);if n==6174||n==0{break;}}seq
}

/// Sudoku solver (backtracking).
pub fn solve_sudoku(grid:&mut [[u8;9];9])->bool{
    for i in 0..9{for j in 0..9{if grid[i][j]==0{
        for num in 1..=9{if is_valid_placement(grid,i,j,num){grid[i][j]=num;
            if solve_sudoku(grid){return true;}grid[i][j]=0;}}return false;}}}true
}

fn is_valid_placement(grid:&[[u8;9];9],row:usize,col:usize,num:u8)->bool{
    for k in 0..9{if grid[row][k]==num||grid[k][col]==num{return false;}}
    let(br,bc)=(row/3*3,col/3*3);
    for i in 0..3{for j in 0..3{if grid[br+i][bc+j]==num{return false;}}}true
}

/// N-Queens solver: returns all solutions.
pub fn n_queens(n:usize)->Vec<Vec<usize>>{
    let mut solutions=Vec::new();let mut board=vec![0;n];
    queens_backtrack(&mut board,0,n,&mut solutions);solutions
}

fn queens_backtrack(board:&mut Vec<usize>,row:usize,n:usize,solutions:&mut Vec<Vec<usize>>){
    if row==n{solutions.push(board.clone());return;}
    for col in 0..n{if queens_safe(board,row,col){board[row]=col;queens_backtrack(board,row+1,n,solutions);}}
}

fn queens_safe(board:&[usize],row:usize,col:usize)->bool{
    for i in 0..row{if board[i]==col||board[i] as i32-col as i32==(i as i32-row as i32)||board[i] as i32-col as i32==(row as i32-i as i32){return false;}}true
}

/// Tower of Hanoi: returns moves as (from, to) pairs.
pub fn hanoi(n:u32)->(Vec<(u8,u8)>,u32){
    let mut moves=Vec::new();fn solve(n:u32,from:u8,to:u8,aux:u8,moves:&mut Vec<(u8,u8)>){
        if n==0{return;}solve(n-1,from,aux,to,moves);moves.push((from,to));solve(n-1,aux,to,from,moves);}
    solve(n,0,2,1,&mut moves);let count=moves.len()as u32;(moves,count)
}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_magic(){let sq=magic_square(3);assert!(is_magic_square(&sq));}
    #[test]fn test_collatz(){assert_eq!(*collatz(6).last().unwrap(),1);}
    #[test]fn test_happy(){assert!(is_happy(7));assert!(!is_happy(4));}
    #[test]fn test_narcissistic(){assert!(is_narcissistic(153));}
    #[test]fn test_queens(){assert_eq!(n_queens(8).len(),92);}
    #[test]fn test_hanoi(){let(moves,count)=hanoi(3);assert_eq!(count,7);}
}
