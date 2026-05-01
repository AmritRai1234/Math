//! Network flow: max flow, min cut, Ford-Fulkerson, bipartite matching.

use std::collections::VecDeque;

/// Max flow using Ford-Fulkerson (BFS/Edmonds-Karp).
pub fn max_flow(capacity:&[Vec<f64>],source:usize,sink:usize)->f64{
    let n=capacity.len();let mut residual:Vec<Vec<f64>>=capacity.to_vec();
    let mut total_flow=0.0;
    loop{// BFS to find augmenting path
        let mut parent=vec![None::<usize>;n];let mut visited=vec![false;n];
        let mut queue=VecDeque::new();queue.push_back(source);visited[source]=true;
        while let Some(u)=queue.pop_front(){for v in 0..n{if !visited[v]&&residual[u][v]>1e-10{
            visited[v]=true;parent[v]=Some(u);queue.push_back(v);if v==sink{break;}}}}
        if !visited[sink]{break;} // no more augmenting paths
        // Find bottleneck
        let mut path_flow=f64::INFINITY;let mut v=sink;
        while v!=source{let u=parent[v].unwrap();path_flow=path_flow.min(residual[u][v]);v=u;}
        // Update residual
        v=sink;while v!=source{let u=parent[v].unwrap();residual[u][v]-=path_flow;residual[v][u]+=path_flow;v=u;}
        total_flow+=path_flow;}
    total_flow
}

/// Min cut: returns the set of vertices on the source side after max flow.
pub fn min_cut(capacity:&[Vec<f64>],source:usize,sink:usize)->Vec<usize>{
    let n=capacity.len();let mut residual:Vec<Vec<f64>>=capacity.to_vec();
    loop{let mut parent=vec![None::<usize>;n];let mut visited=vec![false;n];
        let mut queue=VecDeque::new();queue.push_back(source);visited[source]=true;
        while let Some(u)=queue.pop_front(){for v in 0..n{if !visited[v]&&residual[u][v]>1e-10{
            visited[v]=true;parent[v]=Some(u);queue.push_back(v);}}}
        if !visited[sink]{// BFS one more time on final residual to find reachable
            let mut reachable=vec![false;n];let mut q=VecDeque::new();q.push_back(source);reachable[source]=true;
            while let Some(u)=q.pop_front(){for v in 0..n{if !reachable[v]&&residual[u][v]>1e-10{reachable[v]=true;q.push_back(v);}}}
            return(0..n).filter(|&i|reachable[i]).collect();}
        let mut path_flow=f64::INFINITY;let mut v=sink;
        while v!=source{let u=parent[v].unwrap();path_flow=path_flow.min(residual[u][v]);v=u;}
        v=sink;while v!=source{let u=parent[v].unwrap();residual[u][v]-=path_flow;residual[v][u]+=path_flow;v=u;}}
}

/// Bipartite matching using max flow. Returns number of matched pairs.
pub fn bipartite_matching(adj:&[Vec<usize>],n_left:usize,n_right:usize)->usize{
    let n=n_left+n_right+2;let source=n-2;let sink=n-1;
    let mut cap=vec![vec![0.0;n];n];
    for i in 0..n_left{cap[source][i]=1.0;for &j in &adj[i]{cap[i][n_left+j]=1.0;}}
    for j in 0..n_right{cap[n_left+j][sink]=1.0;}
    max_flow(&cap,source,sink)as usize
}

/// Multi-commodity flow feasibility check (simplified: just sum capacities).
pub fn total_capacity(capacity:&[Vec<f64>],source:usize)->f64{capacity[source].iter().sum()}

/// Shortest augmenting path (number of edges via BFS).
pub fn shortest_path_length(capacity:&[Vec<f64>],source:usize,sink:usize)->Option<usize>{
    let n=capacity.len();let mut dist=vec![usize::MAX;n];let mut queue=VecDeque::new();
    dist[source]=0;queue.push_back(source);
    while let Some(u)=queue.pop_front(){for v in 0..n{if dist[v]==usize::MAX&&capacity[u][v]>1e-10{
        dist[v]=dist[u]+1;queue.push_back(v);if v==sink{return Some(dist[v]);}}}}None
}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_max_flow(){// Simple graph: s→a→t with cap 3, s→b→t with cap 2
        let mut cap=vec![vec![0.0;4];4]; // 0=s, 1=a, 2=b, 3=t
        cap[0][1]=3.0;cap[0][2]=2.0;cap[1][3]=3.0;cap[2][3]=2.0;
        assert!((max_flow(&cap,0,3)-5.0).abs()<1e-10);}
    #[test]fn test_bipartite(){let adj=vec![vec![0,1],vec![0],vec![1,2]];
        assert_eq!(bipartite_matching(&adj,3,3),3);}
    #[test]fn test_min_cut(){let mut cap=vec![vec![0.0;4];4];cap[0][1]=3.0;cap[0][2]=2.0;cap[1][3]=3.0;cap[2][3]=2.0;
        let cut=min_cut(&cap,0,3);assert!(cut.contains(&0));}
}
