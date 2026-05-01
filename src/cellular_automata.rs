//! Cellular automata: elementary automata, Game of Life.

/// Elementary cellular automaton (1D, rule 0-255).
pub fn elementary_automaton(rule:u8,width:usize,steps:usize)->Vec<Vec<u8>>{
    let mut grid=vec![vec![0u8;width]];grid[0][width/2]=1;
    for _ in 0..steps{let prev=grid.last().unwrap();let mut row=vec![0u8;width];
        for i in 0..width{let left=if i>0{prev[i-1]}else{0};let center=prev[i];let right=if i+1<width{prev[i+1]}else{0};
            let pattern=(left<<2)|(center<<1)|right;row[i]=(rule>>pattern)&1;}
        grid.push(row);}grid
}

/// Game of Life step on a 2D grid.
pub fn life_step(grid:&[Vec<u8>])->Vec<Vec<u8>>{
    let rows=grid.len();let cols=grid[0].len();
    let mut next=vec![vec![0u8;cols];rows];
    for i in 0..rows{for j in 0..cols{
        let mut neighbors=0u8;
        for di in -1i32..=1{for dj in -1i32..=1{if di==0&&dj==0{continue;}
            let ni=(i as i32+di).rem_euclid(rows as i32)as usize;
            let nj=(j as i32+dj).rem_euclid(cols as i32)as usize;
            neighbors+=grid[ni][nj];}}
        next[i][j]=match(grid[i][j],neighbors){(1,2)|(1,3)=>1,(0,3)=>1,_=>0};}}next
}

/// Run Game of Life for n steps.
pub fn life_run(initial:&[Vec<u8>],steps:usize)->Vec<Vec<Vec<u8>>>{
    let mut history=vec![initial.to_vec()];let mut grid=initial.to_vec();
    for _ in 0..steps{grid=life_step(&grid);history.push(grid.clone());}history
}

/// Count live cells.
pub fn life_population(grid:&[Vec<u8>])->usize{
    grid.iter().flat_map(|row|row.iter()).filter(|&&c|c==1).count()
}

/// Create glider pattern.
pub fn glider(rows:usize,cols:usize)->Vec<Vec<u8>>{
    let mut grid=vec![vec![0u8;cols];rows];
    grid[0][1]=1;grid[1][2]=1;grid[2][0]=1;grid[2][1]=1;grid[2][2]=1;grid
}

/// Create blinker pattern.
pub fn blinker(rows:usize,cols:usize)->Vec<Vec<u8>>{
    let mut grid=vec![vec![0u8;cols];rows];let r=rows/2;let c=cols/2;
    grid[r-1][c]=1;grid[r][c]=1;grid[r+1][c]=1;grid
}

/// Langton's Ant: simulate for n steps. Returns (grid, ant_row, ant_col, ant_dir).
pub fn langtons_ant(size:usize,steps:usize)->(Vec<Vec<u8>>,usize,usize,u8){
    let mut grid=vec![vec![0u8;size];size];
    let(mut r,mut c,mut dir)=(size/2,size/2,0u8);// 0=up,1=right,2=down,3=left
    for _ in 0..steps{if grid[r][c]==0{dir=(dir+1)%4;grid[r][c]=1;}else{dir=(dir+3)%4;grid[r][c]=0;}
        match dir{0=>r=r.wrapping_sub(1)%size,1=>c=(c+1)%size,2=>r=(r+1)%size,3=>c=c.wrapping_sub(1)%size,_=>{}};}
    (grid,r,c,dir)
}

/// Wire World automaton step.
pub fn wireworld_step(grid:&[Vec<u8>])->Vec<Vec<u8>>{
    let rows=grid.len();let cols=grid[0].len();
    let mut next=vec![vec![0u8;cols];rows];
    for i in 0..rows{for j in 0..cols{next[i][j]=match grid[i][j]{
        0=>0, // empty
        1=>2, // electron head → tail
        2=>3, // tail → conductor
        3=>{  // conductor → head if 1 or 2 neighboring heads
            let mut heads=0;for di in -1i32..=1{for dj in -1i32..=1{if di==0&&dj==0{continue;}
                let ni=i as i32+di;let nj=j as i32+dj;
                if ni>=0&&ni<rows as i32&&nj>=0&&nj<cols as i32&&grid[ni as usize][nj as usize]==1{heads+=1;}}}
            if heads==1||heads==2{1}else{3}},_=>0};}}next
}

/// Rule 110 (Turing complete elementary automaton).
pub fn rule110(width:usize,steps:usize)->Vec<Vec<u8>>{elementary_automaton(110,width,steps)}

/// Rule 30 (chaotic elementary automaton).
pub fn rule30(width:usize,steps:usize)->Vec<Vec<u8>>{elementary_automaton(30,width,steps)}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_elementary(){let g=elementary_automaton(30,11,5);assert_eq!(g.len(),6);assert_eq!(g[0][5],1);}
    #[test]fn test_blinker(){let g=blinker(5,5);let g2=life_step(&g);let g3=life_step(&g2);
        assert_eq!(g,g3);} // blinker has period 2
    #[test]fn test_population(){let g=glider(10,10);assert_eq!(life_population(&g),5);}
    #[test]fn test_langton(){let(grid,_,_,_)=langtons_ant(21,100);let pop:usize=grid.iter().flat_map(|r|r.iter()).filter(|&&c|c==1).count();assert!(pop>0);}
}
