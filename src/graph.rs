//! Graph theory: adjacency, BFS, DFS, Dijkstra, Kruskal, topological sort.

use std::collections::{BinaryHeap,HashMap,HashSet,VecDeque};
use std::cmp::Reverse;

/// Weighted directed graph using adjacency list.
#[derive(Clone,Debug)]
pub struct Graph{pub adj:Vec<Vec<(usize,f64)>>,pub n:usize}

impl Graph{
    pub fn new(n:usize)->Self{Self{adj:vec![vec![];n],n}}
    pub fn add_edge(&mut self,u:usize,v:usize,w:f64){self.adj[u].push((v,w));}
    pub fn add_undirected(&mut self,u:usize,v:usize,w:f64){self.adj[u].push((v,w));self.adj[v].push((u,w));}

    /// BFS from source. Returns parent array (-1 = unvisited).
    pub fn bfs(&self,src:usize)->Vec<i64>{
        let mut parent=vec![-1i64;self.n];parent[src]=-2;
        let mut q=VecDeque::new();q.push_back(src);
        while let Some(u)=q.pop_front(){for &(v,_) in &self.adj[u]{if parent[v]==-1{parent[v]=u as i64;q.push_back(v);}}}
        parent
    }

    /// BFS shortest path (unweighted).
    pub fn bfs_path(&self,src:usize,dst:usize)->Option<Vec<usize>>{
        let parent=self.bfs(src);if parent[dst]==-1{return None;}
        let mut path=vec![dst];let mut cur=dst;
        while cur!=src{cur=parent[cur]as usize;path.push(cur);}
        path.reverse();Some(path)
    }

    /// DFS from source. Returns visit order.
    pub fn dfs(&self,src:usize)->Vec<usize>{
        let mut visited=vec![false;self.n];let mut order=Vec::new();
        self.dfs_helper(src,&mut visited,&mut order);order
    }
    fn dfs_helper(&self,u:usize,visited:&mut Vec<bool>,order:&mut Vec<usize>){
        visited[u]=true;order.push(u);
        for &(v,_) in &self.adj[u]{if !visited[v]{self.dfs_helper(v,visited,order);}}
    }

    /// Dijkstra's shortest path. Returns (distances, parents).
    pub fn dijkstra(&self,src:usize)->(Vec<f64>,Vec<i64>){
        let mut dist=vec![f64::INFINITY;self.n];let mut parent=vec![-1i64;self.n];
        dist[src]=0.0;
        let mut heap=BinaryHeap::new();heap.push(Reverse((OrdF64(0.0),src)));
        while let Some(Reverse((OrdF64(d),u)))=heap.pop(){
            if d>dist[u]{continue;}
            for &(v,w) in &self.adj[u]{let nd=d+w;if nd<dist[v]{dist[v]=nd;parent[v]=u as i64;heap.push(Reverse((OrdF64(nd),v)));}}
        }
        (dist,parent)
    }

    /// Reconstruct path from Dijkstra parents.
    pub fn dijkstra_path(&self,src:usize,dst:usize)->Option<(f64,Vec<usize>)>{
        let(dist,parent)=self.dijkstra(src);
        if dist[dst]==f64::INFINITY{return None;}
        let mut path=vec![dst];let mut cur=dst;
        while cur!=src{cur=parent[cur]as usize;path.push(cur);}
        path.reverse();Some((dist[dst],path))
    }

    /// Topological sort (Kahn's algorithm). Returns None if cycle exists.
    pub fn topological_sort(&self)->Option<Vec<usize>>{
        let mut indeg=vec![0usize;self.n];
        for u in 0..self.n{for &(v,_) in &self.adj[u]{indeg[v]+=1;}}
        let mut q:VecDeque<usize>=indeg.iter().enumerate().filter(|(_,&d)|d==0).map(|(i,_)|i).collect();
        let mut order=Vec::new();
        while let Some(u)=q.pop_front(){order.push(u);for &(v,_) in &self.adj[u]{indeg[v]-=1;if indeg[v]==0{q.push_back(v);}}}
        if order.len()==self.n{Some(order)}else{None}
    }

    /// Detect cycle in directed graph.
    pub fn has_cycle(&self)->bool{self.topological_sort().is_none()}

    /// Connected components (undirected). Returns component ID per node.
    pub fn connected_components(&self)->Vec<usize>{
        let mut comp=vec![usize::MAX;self.n];let mut id=0;
        for s in 0..self.n{if comp[s]==usize::MAX{
            let mut q=VecDeque::new();q.push_back(s);comp[s]=id;
            while let Some(u)=q.pop_front(){for &(v,_) in &self.adj[u]{if comp[v]==usize::MAX{comp[v]=id;q.push_back(v);}}}
            id+=1;}}
        comp
    }

    /// Number of connected components.
    pub fn num_components(&self)->usize{
        let comp=self.connected_components();*comp.iter().max().unwrap_or(&0)+1
    }

    /// Kruskal's MST. Returns edges (u, v, weight) and total weight.
    pub fn kruskal_mst(&self)->(Vec<(usize,usize,f64)>,f64){
        let mut edges=Vec::new();
        for u in 0..self.n{for &(v,w) in &self.adj[u]{if u<v{edges.push((u,v,w));}}}
        edges.sort_by(|a,b|a.2.partial_cmp(&b.2).unwrap());
        let mut uf=UnionFind::new(self.n);let mut mst=Vec::new();let mut total=0.0;
        for(u,v,w) in edges{if uf.find(u)!=uf.find(v){uf.union(u,v);mst.push((u,v,w));total+=w;}}
        (mst,total)
    }

    /// In-degree of each node.
    pub fn in_degrees(&self)->Vec<usize>{
        let mut deg=vec![0;self.n];for u in 0..self.n{for &(v,_) in &self.adj[u]{deg[v]+=1;}}deg
    }

    /// Out-degree of each node.
    pub fn out_degrees(&self)->Vec<usize>{self.adj.iter().map(|e|e.len()).collect()}
}

#[derive(Clone,Copy,PartialEq,PartialOrd)]
struct OrdF64(f64);
impl Eq for OrdF64{}
impl Ord for OrdF64{fn cmp(&self,o:&Self)->std::cmp::Ordering{self.0.partial_cmp(&o.0).unwrap_or(std::cmp::Ordering::Equal)}}

struct UnionFind{parent:Vec<usize>,rank:Vec<usize>}
impl UnionFind{
    fn new(n:usize)->Self{Self{parent:(0..n).collect(),rank:vec![0;n]}}
    fn find(&mut self,x:usize)->usize{if self.parent[x]!=x{self.parent[x]=self.find(self.parent[x]);}self.parent[x]}
    fn union(&mut self,x:usize,y:usize){let(rx,ry)=(self.find(x),self.find(y));if rx==ry{return;}
        if self.rank[rx]<self.rank[ry]{self.parent[rx]=ry;}else{self.parent[ry]=rx;if self.rank[rx]==self.rank[ry]{self.rank[rx]+=1;}}}
}

/// Create adjacency matrix from graph.
pub fn to_adjacency_matrix(g:&Graph)->Vec<Vec<f64>>{
    let mut m=vec![vec![0.0;g.n];g.n];
    for u in 0..g.n{for &(v,w) in &g.adj[u]{m[u][v]=w;}}m
}

/// Floyd-Warshall all-pairs shortest paths.
pub fn floyd_warshall(n:usize,edges:&[(usize,usize,f64)])->Vec<Vec<f64>>{
    let mut dist=vec![vec![f64::INFINITY;n];n];
    for i in 0..n{dist[i][i]=0.0;}
    for &(u,v,w) in edges{dist[u][v]=w;}
    for k in 0..n{for i in 0..n{for j in 0..n{let d=dist[i][k]+dist[k][j];if d<dist[i][j]{dist[i][j]=d;}}}}
    dist
}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_bfs(){let mut g=Graph::new(4);g.add_undirected(0,1,1.0);g.add_undirected(1,2,1.0);g.add_undirected(2,3,1.0);
        let path=g.bfs_path(0,3).unwrap();assert_eq!(path,vec![0,1,2,3]);}
    #[test]fn test_dijkstra(){let mut g=Graph::new(3);g.add_edge(0,1,1.0);g.add_edge(1,2,2.0);g.add_edge(0,2,10.0);
        let(d,_)=g.dijkstra(0);assert!((d[2]-3.0).abs()<1e-10);}
    #[test]fn test_topo(){let mut g=Graph::new(4);g.add_edge(0,1,1.0);g.add_edge(0,2,1.0);g.add_edge(1,3,1.0);g.add_edge(2,3,1.0);
        let order=g.topological_sort().unwrap();assert_eq!(order[0],0);assert_eq!(*order.last().unwrap(),3);}
    #[test]fn test_kruskal(){let mut g=Graph::new(4);g.add_undirected(0,1,1.0);g.add_undirected(1,2,2.0);g.add_undirected(2,3,3.0);g.add_undirected(0,3,4.0);
        let(_,total)=g.kruskal_mst();assert!((total-6.0).abs()<1e-10);}
    #[test]fn test_cycle(){let mut g=Graph::new(3);g.add_edge(0,1,1.0);g.add_edge(1,2,1.0);g.add_edge(2,0,1.0);assert!(g.has_cycle());}
}
