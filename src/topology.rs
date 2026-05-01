//! Computational topology: simplicial complexes, Euler characteristic, Betti numbers.

/// Simplicial complex represented as sets of simplices (vertices, edges, faces, etc.).
#[derive(Clone,Debug)]
pub struct SimplicialComplex{pub simplices:Vec<Vec<Vec<usize>>>,pub dim:usize}

impl SimplicialComplex{
    /// Create from maximal simplices.
    pub fn from_facets(facets:&[Vec<usize>])->Self{
        let max_dim=facets.iter().map(|f|f.len()).max().unwrap_or(0);
        let mut simplices=vec![Vec::new();max_dim];
        let mut seen=std::collections::HashSet::new();
        for facet in facets{Self::add_faces(facet,&mut simplices,&mut seen);}
        Self{dim:max_dim.saturating_sub(1),simplices}
    }

    fn add_faces(simplex:&[usize],simplices:&mut Vec<Vec<Vec<usize>>>,seen:&mut std::collections::HashSet<Vec<usize>>){
        let n=simplex.len();
        // Generate all subsets (faces)
        for mask in 1..(1u64<<n){let mut face:Vec<usize>=Vec::new();
            for i in 0..n{if mask&(1<<i)!=0{face.push(simplex[i]);}}
            face.sort();if seen.insert(face.clone()){let dim=face.len()-1;
                if dim<simplices.len(){simplices[dim].push(face);}}}
    }

    /// Number of k-simplices.
    pub fn count(&self,k:usize)->usize{if k<self.simplices.len(){self.simplices[k].len()}else{0}}

    /// Euler characteristic χ = Σ(-1)^k * f_k.
    pub fn euler_characteristic(&self)->i64{
        self.simplices.iter().enumerate().map(|(k,s)|if k%2==0{s.len()as i64}else{-(s.len()as i64)}).sum()
    }

    /// f-vector: (vertices, edges, faces, ...).
    pub fn f_vector(&self)->Vec<usize>{self.simplices.iter().map(|s|s.len()).collect()}
}

/// Euler characteristic for polyhedra: V - E + F.
pub fn euler_polyhedron(vertices:usize,edges:usize,faces:usize)->i64{vertices as i64-edges as i64+faces as i64}

/// Euler characteristic for surface with g handles: 2 - 2g.
pub fn euler_surface(genus:usize)->i64{2-2*genus as i64}

/// Check Euler formula for convex polyhedron (should be 2).
pub fn is_convex_polyhedron(v:usize,e:usize,f:usize)->bool{euler_polyhedron(v,e,f)==2}

/// Genus from Euler characteristic: g = (2 - χ) / 2.
pub fn genus_from_euler(chi:i64)->f64{(2-chi)as f64/2.0}

/// Graph planarity check (simplified): K5 and K3,3 are not planar.
/// Uses Euler's formula: for planar graph, E ≤ 3V - 6.
pub fn is_planar_bound(vertices:usize,edges:usize)->bool{
    if vertices<3{return true;}edges<=3*vertices-6
}

/// Compute boundary matrix for chain complex (from k-simplices to (k-1)-simplices).
pub fn boundary_matrix(k_simplices:&[Vec<usize>],k_minus_1_simplices:&[Vec<usize>])->Vec<Vec<i8>>{
    let m=k_minus_1_simplices.len();let n=k_simplices.len();
    let mut matrix=vec![vec![0i8;n];m];
    for(j,simplex) in k_simplices.iter().enumerate(){
        for(i,face) in simplex.iter().enumerate(){
            let mut boundary=simplex.clone();boundary.remove(i);boundary.sort();
            if let Some(pos)=k_minus_1_simplices.iter().position(|s|*s==boundary){
                matrix[pos][j]=if i%2==0{1}else{-1};}}}
    matrix
}

/// Compute Betti number β₀ (connected components) from vertex and edge lists.
pub fn betti_0(n_vertices:usize,edges:&[(usize,usize)])->usize{
    let mut parent:Vec<usize>=(0..n_vertices).collect();
    fn find(parent:&mut Vec<usize>,x:usize)->usize{if parent[x]!=x{parent[x]=find(parent,parent[x]);}parent[x]}
    for &(u,v) in edges{let(pu,pv)=(find(&mut parent,u),find(&mut parent,v));if pu!=pv{parent[pu]=pv;}}
    let mut components=std::collections::HashSet::new();
    for i in 0..n_vertices{components.insert(find(&mut parent,i));}components.len()
}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_tetrahedron(){let sc=SimplicialComplex::from_facets(&[vec![0,1,2],vec![0,1,3],vec![0,2,3],vec![1,2,3]]);
        assert_eq!(sc.euler_characteristic(),2);} // sphere has χ=2
    #[test]fn test_euler_cube(){assert_eq!(euler_polyhedron(8,12,6),2);}
    #[test]fn test_euler_surface(){assert_eq!(euler_surface(0),2);assert_eq!(euler_surface(1),0);}
    #[test]fn test_betti0(){assert_eq!(betti_0(4,&[(0,1),(2,3)]),2);} // two components
    #[test]fn test_planar(){assert!(is_planar_bound(4,6));assert!(!is_planar_bound(5,11));}
}
