//! Clustering algorithms: k-means, DBSCAN, hierarchical.

use crate::powers;

/// K-means clustering. Returns (assignments, centroids).
pub fn kmeans(data:&[Vec<f64>],k:usize,max_iter:usize)->(Vec<usize>,Vec<Vec<f64>>){
    let n=data.len();let d=data[0].len();
    // Initialize centroids from first k points
    let mut centroids:Vec<Vec<f64>>=data.iter().take(k).cloned().collect();
    let mut assign=vec![0usize;n];
    for _ in 0..max_iter{
        // Assign
        let mut changed=false;
        for i in 0..n{let mut best=0;let mut best_d=f64::INFINITY;
            for c in 0..k{let dist=euclidean_dist(&data[i],&centroids[c]);if dist<best_d{best_d=dist;best=c;}}
            if assign[i]!=best{changed=true;assign[i]=best;}}
        if !changed{break;}
        // Update centroids
        let mut sums=vec![vec![0.0;d];k];let mut counts=vec![0usize;k];
        for i in 0..n{let c=assign[i];counts[c]+=1;for j in 0..d{sums[c][j]+=data[i][j];}}
        for c in 0..k{if counts[c]>0{for j in 0..d{centroids[c][j]=sums[c][j]/counts[c]as f64;}}}
    }
    (assign,centroids)
}

/// K-means++ initialization.
pub fn kmeans_pp(data:&[Vec<f64>],k:usize,max_iter:usize)->(Vec<usize>,Vec<Vec<f64>>){
    let n=data.len();let mut centroids=vec![data[0].clone()];let mut rng=42u64;
    while centroids.len()<k{let mut dists:Vec<f64>=data.iter().map(|p|{
        centroids.iter().map(|c|euclidean_dist(p,c)).fold(f64::INFINITY,f64::min)}).collect();
        let total:f64=dists.iter().map(|d|d*d).sum();
        rng=crate::noise::xorshift(rng);let target=(rng%10000)as f64/10000.0*total;
        let mut cum=0.0;for i in 0..n{cum+=dists[i]*dists[i];if cum>=target{centroids.push(data[i].clone());break;}}}
    let mut result_data:Vec<Vec<f64>>=data.to_vec();
    // Run regular kmeans with better init
    let mut assign=vec![0usize;n];
    for _ in 0..max_iter{let mut changed=false;
        for i in 0..n{let mut best=0;let mut best_d=f64::INFINITY;
            for c in 0..k{let dist=euclidean_dist(&data[i],&centroids[c]);if dist<best_d{best_d=dist;best=c;}}
            if assign[i]!=best{changed=true;assign[i]=best;}}
        if !changed{break;}
        let d=data[0].len();let mut sums=vec![vec![0.0;d];k];let mut counts=vec![0usize;k];
        for i in 0..n{let c=assign[i];counts[c]+=1;for j in 0..d{sums[c][j]+=data[i][j];}}
        for c in 0..k{if counts[c]>0{for j in 0..d{centroids[c][j]=sums[c][j]/counts[c]as f64;}}}}
    (assign,centroids)
}

/// DBSCAN clustering. Returns cluster labels (-1 = noise).
pub fn dbscan(data:&[Vec<f64>],eps:f64,min_pts:usize)->Vec<i32>{
    let n=data.len();let mut labels=vec![-1i32;n];let mut cluster=0i32;
    for i in 0..n{if labels[i]!=-1{continue;}
        let neighbors=range_query(data,i,eps);
        if neighbors.len()<min_pts{continue;} // noise (might be claimed later)
        labels[i]=cluster;let mut queue=neighbors;let mut qi=0;
        while qi<queue.len(){let j=queue[qi];qi+=1;
            if labels[j]==-1{labels[j]=cluster;}
            else if labels[j]!=-1{continue;} // skip if already in a cluster
            labels[j]=cluster;let nb=range_query(data,j,eps);
            if nb.len()>=min_pts{for &k in &nb{if !queue.contains(&k){queue.push(k);}}}}
        cluster+=1;}
    labels
}

fn range_query(data:&[Vec<f64>],idx:usize,eps:f64)->Vec<usize>{
    (0..data.len()).filter(|&i|euclidean_dist(&data[idx],&data[i])<=eps).collect()
}

/// Hierarchical agglomerative clustering (single linkage). Returns merge history.
pub fn hierarchical(data:&[Vec<f64>],k:usize)->Vec<usize>{
    let n=data.len();let mut assign:Vec<usize>=(0..n).collect();let mut num_clusters=n;
    while num_clusters>k{let mut best_i=0;let mut best_j=1;let mut best_d=f64::INFINITY;
        for i in 0..n{for j in i+1..n{if assign[i]!=assign[j]{
            let d=euclidean_dist(&data[i],&data[j]);if d<best_d{best_d=d;best_i=i;best_j=j;}}}}
        let old=assign[best_j];let new=assign[best_i];
        for i in 0..n{if assign[i]==old{assign[i]=new;}}num_clusters-=1;}
    // Renumber clusters 0..k-1
    let mut seen=std::collections::HashMap::new();let mut id=0;
    for a in &mut assign{let entry=seen.entry(*a).or_insert_with(||{let v=id;id+=1;v});*a=*entry;}
    assign
}

/// Silhouette score for clustering quality (-1 to 1, higher is better).
pub fn silhouette_score(data:&[Vec<f64>],labels:&[usize])->f64{
    let n=data.len();if n<2{return 0.0;}
    let mut total=0.0;
    for i in 0..n{let ci=labels[i];
        let mut a_sum=0.0;let mut a_count=0;
        for j in 0..n{if j!=i&&labels[j]==ci{a_sum+=euclidean_dist(&data[i],&data[j]);a_count+=1;}}
        let a=if a_count>0{a_sum/a_count as f64}else{0.0};
        let mut b=f64::INFINITY;
        let k=*labels.iter().max().unwrap()+1;
        for c in 0..k{if c!=ci{let mut sum=0.0;let mut count=0;
            for j in 0..n{if labels[j]==c{sum+=euclidean_dist(&data[i],&data[j]);count+=1;}}
            if count>0{let avg=sum/count as f64;if avg<b{b=avg;}}}}
        let s=if a.max(b)>0.0{(b-a)/a.max(b)}else{0.0};total+=s;}
    total/n as f64
}

fn euclidean_dist(a:&[f64],b:&[f64])->f64{
    powers::sqrt(a.iter().zip(b.iter()).map(|(&x,&y)|(x-y)*(x-y)).sum())
}

/// Inertia (within-cluster sum of squares).
pub fn inertia(data:&[Vec<f64>],labels:&[usize],centroids:&[Vec<f64>])->f64{
    data.iter().zip(labels.iter()).map(|(p,&l)|{let c=&centroids[l];p.iter().zip(c.iter()).map(|(&a,&b)|(a-b)*(a-b)).sum::<f64>()}).sum()
}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_kmeans(){let data=vec![vec![0.0,0.0],vec![0.1,0.1],vec![10.0,10.0],vec![10.1,10.1]];
        let(assign,_)=kmeans(&data,2,100);assert_eq!(assign[0],assign[1]);assert_eq!(assign[2],assign[3]);assert_ne!(assign[0],assign[2]);}
    #[test]fn test_dbscan(){let data=vec![vec![0.0],vec![0.1],vec![0.2],vec![10.0],vec![10.1]];
        let labels=dbscan(&data,0.5,2);assert_eq!(labels[0],labels[1]);assert_eq!(labels[3],labels[4]);}
}
