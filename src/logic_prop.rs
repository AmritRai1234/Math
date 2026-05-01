//! Symbolic logic: truth tables, propositional calculus, SAT solver.

/// Propositional formula.
#[derive(Clone,Debug)]
pub enum Prop{Var(usize),Not(Box<Prop>),And(Box<Prop>,Box<Prop>),Or(Box<Prop>,Box<Prop>),
    Implies(Box<Prop>,Box<Prop>),Iff(Box<Prop>,Box<Prop>),True,False}

impl Prop{
    pub fn var(i:usize)->Self{Prop::Var(i)}
    pub fn not(p:Prop)->Self{Prop::Not(Box::new(p))}
    pub fn and(a:Prop,b:Prop)->Self{Prop::And(Box::new(a),Box::new(b))}
    pub fn or(a:Prop,b:Prop)->Self{Prop::Or(Box::new(a),Box::new(b))}
    pub fn implies(a:Prop,b:Prop)->Self{Prop::Implies(Box::new(a),Box::new(b))}
    pub fn iff(a:Prop,b:Prop)->Self{Prop::Iff(Box::new(a),Box::new(b))}

    /// Evaluate formula under assignment.
    pub fn eval(&self,vals:&[bool])->bool{match self{
        Prop::Var(i)=>vals[*i],Prop::Not(p)=>!p.eval(vals),
        Prop::And(a,b)=>a.eval(vals)&&b.eval(vals),Prop::Or(a,b)=>a.eval(vals)||b.eval(vals),
        Prop::Implies(a,b)=>!a.eval(vals)||b.eval(vals),Prop::Iff(a,b)=>a.eval(vals)==b.eval(vals),
        Prop::True=>true,Prop::False=>false}}

    /// Count variables.
    pub fn num_vars(&self)->usize{let mut max=0;self.collect_vars(&mut max);max+1}
    fn collect_vars(&self,max:&mut usize){match self{
        Prop::Var(i)=>{if *i>*max{*max=*i;}},Prop::Not(p)=>p.collect_vars(max),
        Prop::And(a,b)|Prop::Or(a,b)|Prop::Implies(a,b)|Prop::Iff(a,b)=>{a.collect_vars(max);b.collect_vars(max);},_=>{}}}
}

/// Generate truth table: returns all assignments and their results.
pub fn truth_table(formula:&Prop,n_vars:usize)->Vec<(Vec<bool>,bool)>{
    let mut results=Vec::new();
    for mask in 0..(1u64<<n_vars){let vals:Vec<bool>=(0..n_vars).map(|i|mask&(1<<i)!=0).collect();
        let result=formula.eval(&vals);results.push((vals,result));}results
}

/// Check if formula is a tautology.
pub fn is_tautology(formula:&Prop)->bool{
    let n=formula.num_vars();truth_table(formula,n).iter().all(|(_,r)|*r)
}

/// Check if formula is satisfiable.
pub fn is_satisfiable(formula:&Prop)->bool{
    let n=formula.num_vars();truth_table(formula,n).iter().any(|(_,r)|*r)
}

/// Check if formula is a contradiction.
pub fn is_contradiction(formula:&Prop)->bool{!is_satisfiable(formula)}

/// Find satisfying assignment.
pub fn find_satisfying(formula:&Prop)->Option<Vec<bool>>{
    let n=formula.num_vars();truth_table(formula,n).into_iter().find(|(_,r)|*r).map(|(v,_)|v)
}

/// Check logical equivalence.
pub fn are_equivalent(a:&Prop,b:&Prop)->bool{is_tautology(&Prop::iff(a.clone(),b.clone()))}

/// CNF clause: list of literals (positive = var, negative = ¬var).
pub type Clause=Vec<i32>;

/// DPLL SAT solver.
pub fn dpll_sat(clauses:&[Clause],n_vars:usize)->Option<Vec<bool>>{
    let mut assignment=vec![None::<bool>;n_vars+1];
    if dpll_recurse(clauses,&mut assignment){Some((1..=n_vars).map(|i|assignment[i].unwrap_or(false)).collect())}else{None}
}

fn dpll_recurse(clauses:&[Clause],assignment:&mut Vec<Option<bool>>)->bool{
    // Unit propagation
    let mut clauses=clauses.to_vec();
    loop{let mut changed=false;
        for i in 0..clauses.len(){if clauses[i].len()==1{let lit=clauses[i][0];
            let var=lit.unsigned_abs()as usize;let val=lit>0;
            assignment[var]=Some(val);clauses=propagate(&clauses,lit);changed=true;break;}}
        if !changed{break;}}
    // Check
    if clauses.is_empty(){return true;}
    if clauses.iter().any(|c|c.is_empty()){return false;}
    // Choose variable
    let var=clauses[0][0].unsigned_abs()as usize;
    // Try true
    assignment[var]=Some(true);
    let pos=propagate(&clauses,var as i32);
    if dpll_recurse(&pos,assignment){return true;}
    // Try false
    assignment[var]=Some(false);
    let neg=propagate(&clauses,-(var as i32));
    if dpll_recurse(&neg,assignment){return true;}
    assignment[var]=None;false
}

fn propagate(clauses:&[Clause],lit:i32)->Vec<Clause>{
    clauses.iter().filter(|c|!c.contains(&lit))
        .map(|c|c.iter().filter(|&&l|l!=-lit).cloned().collect()).collect()
}

/// Count satisfying assignments.
pub fn count_sat(formula:&Prop)->usize{
    let n=formula.num_vars();truth_table(formula,n).iter().filter(|(_,r)|*r).count()
}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_tautology(){let p=Prop::or(Prop::var(0),Prop::not(Prop::var(0)));assert!(is_tautology(&p));}
    #[test]fn test_contradiction(){let p=Prop::and(Prop::var(0),Prop::not(Prop::var(0)));assert!(is_contradiction(&p));}
    #[test]fn test_sat(){let p=Prop::and(Prop::var(0),Prop::var(1));assert!(is_satisfiable(&p));}
    #[test]fn test_dpll(){let clauses=vec![vec![1,2],vec![-1,2],vec![1,-2]];
        let sol=dpll_sat(&clauses,2).unwrap();assert!(sol[0]||sol[1]);}
    #[test]fn test_equiv(){let a=Prop::implies(Prop::var(0),Prop::var(1));
        let b=Prop::or(Prop::not(Prop::var(0)),Prop::var(1));assert!(are_equivalent(&a,&b));}
}
