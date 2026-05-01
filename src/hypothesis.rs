//! Hypothesis testing: t-test, chi-squared, ANOVA, p-values.

use crate::powers;use crate::stats;use crate::arithmetic::abs;

/// One-sample t-test. Returns (t_statistic, degrees_of_freedom).
pub fn t_test_one_sample(data:&[f64],mu0:f64)->(f64,usize){
    let n=data.len();let mean=stats::mean(data);let s=stats::std_dev(data);
    let t=(mean-mu0)/(s/powers::sqrt(n as f64));(t,n-1)
}

/// Two-sample t-test (equal variance assumed). Returns (t_statistic, df).
pub fn t_test_two_sample(a:&[f64],b:&[f64])->(f64,usize){
    let(n1,n2)=(a.len(),b.len());let(m1,m2)=(stats::mean(a),stats::mean(b));
    let s1=stats::variance(a);let s2=stats::variance(b);
    let sp=powers::sqrt(((n1-1)as f64*s1+(n2-1)as f64*s2)/(n1+n2-2)as f64);
    let t=(m1-m2)/(sp*powers::sqrt(1.0/n1 as f64+1.0/n2 as f64));
    (t,n1+n2-2)
}

/// Welch's t-test (unequal variance). Returns (t_statistic, approx_df).
pub fn welch_t_test(a:&[f64],b:&[f64])->(f64,f64){
    let(n1,n2)=(a.len()as f64,b.len()as f64);let(m1,m2)=(stats::mean(a),stats::mean(b));
    let(s1,s2)=(stats::variance(a),stats::variance(b));
    let se=powers::sqrt(s1/n1+s2/n2);let t=(m1-m2)/se;
    let num=(s1/n1+s2/n2)*(s1/n1+s2/n2);
    let den=(s1/n1)*(s1/n1)/(n1-1.0)+(s2/n2)*(s2/n2)/(n2-1.0);
    (t,num/den)
}

/// Paired t-test. Returns (t_statistic, df).
pub fn t_test_paired(a:&[f64],b:&[f64])->(f64,usize){
    let diffs:Vec<f64>=a.iter().zip(b.iter()).map(|(&x,&y)|x-y).collect();
    t_test_one_sample(&diffs,0.0)
}

/// One-way ANOVA. Returns (f_statistic, df_between, df_within).
pub fn anova_one_way(groups:&[Vec<f64>])->(f64,usize,usize){
    let k=groups.len();let n:usize=groups.iter().map(|g|g.len()).sum();
    let grand_mean:f64=groups.iter().flat_map(|g|g.iter()).sum::<f64>()/n as f64;
    let ss_between:f64=groups.iter().map(|g|g.len()as f64*(stats::mean(g)-grand_mean)*(stats::mean(g)-grand_mean)).sum();
    let ss_within:f64=groups.iter().flat_map(|g|g.iter().map(|&x|(x-stats::mean(g))*(x-stats::mean(g)))).sum();
    let df_b=k-1;let df_w=n-k;
    let f=(ss_between/df_b as f64)/(ss_within/df_w as f64);(f,df_b,df_w)
}

/// Chi-squared test for independence. Returns (chi2, df).
pub fn chi_squared_test(observed:&[Vec<f64>])->(f64,usize){
    let rows=observed.len();let cols=observed[0].len();
    let row_sums:Vec<f64>=observed.iter().map(|r|r.iter().sum()).collect();
    let col_sums:Vec<f64>=(0..cols).map(|j|observed.iter().map(|r|r[j]).sum()).collect();
    let total:f64=row_sums.iter().sum();
    let mut chi2=0.0;
    for i in 0..rows{for j in 0..cols{let expected=row_sums[i]*col_sums[j]/total;
        if expected>0.0{chi2+=(observed[i][j]-expected)*(observed[i][j]-expected)/expected;}}}
    (chi2,(rows-1)*(cols-1))
}

/// Chi-squared goodness of fit. Returns (chi2, df).
pub fn chi_squared_gof(observed:&[f64],expected:&[f64])->(f64,usize){
    let chi2:f64=observed.iter().zip(expected.iter()).map(|(&o,&e)|if e>0.0{(o-e)*(o-e)/e}else{0.0}).sum();
    (chi2,observed.len()-1)
}

/// Z-test for proportions. Returns z-statistic.
pub fn z_test_proportion(successes:u64,n:u64,p0:f64)->f64{
    let p_hat=successes as f64/n as f64;let se=powers::sqrt(p0*(1.0-p0)/n as f64);(p_hat-p0)/se
}

/// Mann-Whitney U statistic.
pub fn mann_whitney_u(a:&[f64],b:&[f64])->(f64,f64){
    let(n1,n2)=(a.len(),b.len());
    let u1:f64=a.iter().map(|&x|b.iter().filter(|&&y|x>y).count()as f64+b.iter().filter(|&&y|(x-y).abs()<1e-15).count()as f64*0.5).sum();
    let u2=n1 as f64*n2 as f64-u1;(u1,u2)
}

/// Cohen's d effect size.
pub fn cohens_d(a:&[f64],b:&[f64])->f64{
    let(m1,m2)=(stats::mean(a),stats::mean(b));
    let sp=powers::sqrt((stats::variance(a)+stats::variance(b))/2.0);
    if sp<1e-30{return 0.0;}(m1-m2)/sp
}

/// Pearson's r (correlation coefficient).
pub fn pearson_r(x:&[f64],y:&[f64])->f64{stats::correlation(x,y)}

/// Spearman's rank correlation.
pub fn spearman_rho(x:&[f64],y:&[f64])->f64{
    let rank=|data:&[f64]|->Vec<f64>{let mut indexed:Vec<(usize,f64)>=data.iter().enumerate().map(|(i,&v)|(i,v)).collect();
        indexed.sort_by(|a,b|a.1.partial_cmp(&b.1).unwrap());
        let mut ranks=vec![0.0;data.len()];for(rank,&(i,_)) in indexed.iter().enumerate(){ranks[i]=rank as f64+1.0;}ranks};
    let rx=rank(x);let ry=rank(y);pearson_r(&rx,&ry)
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_t_one(){let data=[5.1,4.9,5.0,5.2,4.8];let(t,df)=t_test_one_sample(&data,5.0);assert_eq!(df,4);}
    #[test]fn test_anova(){let g=vec![vec![1.0,2.0,3.0],vec![4.0,5.0,6.0]];let(f,_,_)=anova_one_way(&g);assert!(f>0.0);}
    #[test]fn test_chi2(){let obs=vec![vec![10.0,20.0],vec![30.0,40.0]];let(chi2,df)=chi_squared_test(&obs);assert_eq!(df,1);}
    #[test]fn test_cohens_d(){let a=[10.0,11.0,12.0];let b=[8.0,9.0,10.0];let d=cohens_d(&a,&b);assert!(d>0.0);}
}
