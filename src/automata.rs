//! Finite automata: DFA, NFA, state machines.

use std::collections::{HashMap,HashSet,VecDeque};

/// Deterministic Finite Automaton.
#[derive(Clone,Debug)]
pub struct DFA{pub states:usize,pub alphabet:Vec<char>,pub transitions:HashMap<(usize,char),usize>,pub start:usize,pub accept:HashSet<usize>}

impl DFA{
    pub fn new(states:usize,alphabet:Vec<char>,start:usize,accept:HashSet<usize>)->Self{
        Self{states,alphabet,transitions:HashMap::new(),start,accept}
    }
    pub fn add_transition(&mut self,from:usize,on:char,to:usize){self.transitions.insert((from,on),to);}

    /// Run DFA on input string. Returns true if accepted.
    pub fn accepts(&self,input:&str)->bool{
        let mut state=self.start;
        for c in input.chars(){match self.transitions.get(&(state,c)){Some(&next)=>state=next,None=>return false,}}
        self.accept.contains(&state)
    }

    /// Minimize DFA using Hopcroft's algorithm (simplified).
    pub fn minimize(&self)->DFA{
        // Remove unreachable states
        let mut reachable=HashSet::new();let mut queue=VecDeque::new();
        reachable.insert(self.start);queue.push_back(self.start);
        while let Some(s)=queue.pop_front(){for &c in &self.alphabet{
            if let Some(&next)=self.transitions.get(&(s,c)){if reachable.insert(next){queue.push_back(next);}}}}
        // Partition into accepting and non-accepting
        let accept:HashSet<usize>=self.accept.iter().filter(|s|reachable.contains(s)).cloned().collect();
        let non_accept:HashSet<usize>=reachable.iter().filter(|s|!accept.contains(s)).cloned().collect();
        let mut partitions:Vec<HashSet<usize>>=vec![];
        if !non_accept.is_empty(){partitions.push(non_accept);}
        if !accept.is_empty(){partitions.push(accept);}
        // Refine partitions
        let mut changed=true;
        while changed{changed=false;let mut new_parts=Vec::new();
            for part in &partitions{let mut splits:HashMap<Vec<usize>,HashSet<usize>>=HashMap::new();
                for &s in part{let key:Vec<usize>=self.alphabet.iter().map(|&c|{
                    let next=self.transitions.get(&(s,c)).cloned().unwrap_or(self.states);
                    partitions.iter().position(|p|p.contains(&next)).unwrap_or(partitions.len())}).collect();
                    splits.entry(key).or_default().insert(s);}
                if splits.len()>1{changed=true;}
                for(_,set)in splits{new_parts.push(set);}}
            partitions=new_parts;}
        // Build minimized DFA
        let state_map:HashMap<usize,usize>=partitions.iter().enumerate()
            .flat_map(|(i,set)|set.iter().map(move|&s|(s,i))).collect();
        let n=partitions.len();
        let new_start=state_map[&self.start];
        let new_accept:HashSet<usize>=partitions.iter().enumerate()
            .filter(|(_,set)|set.iter().any(|s|self.accept.contains(s))).map(|(i,_)|i).collect();
        let mut dfa=DFA::new(n,self.alphabet.clone(),new_start,new_accept);
        for(&(from,c),&to) in &self.transitions{if reachable.contains(&from){
            dfa.add_transition(state_map[&from],c,state_map[&to]);}}
        dfa
    }
}

/// Non-deterministic Finite Automaton.
#[derive(Clone,Debug)]
pub struct NFA{pub states:usize,pub alphabet:Vec<char>,pub transitions:HashMap<(usize,char),HashSet<usize>>,
    pub epsilon:HashMap<usize,HashSet<usize>>,pub start:usize,pub accept:HashSet<usize>}

impl NFA{
    pub fn new(states:usize,alphabet:Vec<char>,start:usize,accept:HashSet<usize>)->Self{
        Self{states,alphabet,transitions:HashMap::new(),epsilon:HashMap::new(),start,accept}
    }
    pub fn add_transition(&mut self,from:usize,on:char,to:usize){self.transitions.entry((from,on)).or_default().insert(to);}
    pub fn add_epsilon(&mut self,from:usize,to:usize){self.epsilon.entry(from).or_default().insert(to);}

    fn epsilon_closure(&self,states:&HashSet<usize>)->HashSet<usize>{
        let mut closure=states.clone();let mut queue:VecDeque<usize>=states.iter().cloned().collect();
        while let Some(s)=queue.pop_front(){if let Some(eps)=self.epsilon.get(&s){
            for &next in eps{if closure.insert(next){queue.push_back(next);}}}}closure
    }

    pub fn accepts(&self,input:&str)->bool{
        let mut current=self.epsilon_closure(&[self.start].iter().cloned().collect());
        for c in input.chars(){let mut next=HashSet::new();
            for &s in &current{if let Some(targets)=self.transitions.get(&(s,c)){next.extend(targets);}}
            current=self.epsilon_closure(&next);if current.is_empty(){return false;}}
        current.iter().any(|s|self.accept.contains(s))
    }

    /// Convert NFA to DFA (subset construction).
    pub fn to_dfa(&self)->DFA{
        let start_set=self.epsilon_closure(&[self.start].iter().cloned().collect());
        let mut dfa_states:Vec<HashSet<usize>>=vec![start_set.clone()];
        let mut queue=VecDeque::new();queue.push_back(0usize);
        let mut dfa=DFA::new(0,self.alphabet.clone(),0,HashSet::new());
        if start_set.iter().any(|s|self.accept.contains(s)){dfa.accept.insert(0);}
        while let Some(idx)=queue.pop_front(){let current=dfa_states[idx].clone();
            for &c in &self.alphabet{let mut next=HashSet::new();
                for &s in &current{if let Some(targets)=self.transitions.get(&(s,c)){next.extend(targets);}}
                let next=self.epsilon_closure(&next);if next.is_empty(){continue;}
                let target_idx=if let Some(pos)=dfa_states.iter().position(|s|*s==next){pos}
                    else{let pos=dfa_states.len();dfa_states.push(next.clone());queue.push_back(pos);
                    if next.iter().any(|s|self.accept.contains(s)){dfa.accept.insert(pos);}pos};
                dfa.add_transition(idx,c,target_idx);}}
        dfa.states=dfa_states.len();dfa
    }
}

/// Mealy machine.
pub struct Mealy{pub transitions:HashMap<(usize,char),(usize,char)>,pub start:usize}
impl Mealy{
    pub fn new(start:usize)->Self{Self{transitions:HashMap::new(),start}}
    pub fn add(&mut self,from:usize,input:char,to:usize,output:char){self.transitions.insert((from,input),(to,output));}
    pub fn run(&self,input:&str)->String{let mut state=self.start;let mut out=String::new();
        for c in input.chars(){if let Some(&(next,o))=self.transitions.get(&(state,c)){state=next;out.push(o);}}out}
}

#[cfg(test)]
mod tests{use super::*;
    #[test]fn test_dfa(){// DFA accepting strings ending in "ab"
        let mut dfa=DFA::new(3,vec!['a','b'],0,[2].iter().cloned().collect());
        dfa.add_transition(0,'a',1);dfa.add_transition(0,'b',0);
        dfa.add_transition(1,'a',1);dfa.add_transition(1,'b',2);
        dfa.add_transition(2,'a',1);dfa.add_transition(2,'b',0);
        assert!(dfa.accepts("ab"));assert!(dfa.accepts("aab"));assert!(!dfa.accepts("ba"));}
    #[test]fn test_nfa(){let mut nfa=NFA::new(3,vec!['a','b'],0,[2].iter().cloned().collect());
        nfa.add_transition(0,'a',0);nfa.add_transition(0,'a',1);nfa.add_transition(1,'b',2);
        assert!(nfa.accepts("ab"));assert!(nfa.accepts("aab"));assert!(!nfa.accepts("b"));}
    #[test]fn test_nfa_to_dfa(){let mut nfa=NFA::new(2,vec!['a','b'],0,[1].iter().cloned().collect());
        nfa.add_transition(0,'a',1);nfa.add_transition(0,'b',0);
        let dfa=nfa.to_dfa();assert!(dfa.accepts("a"));assert!(!dfa.accepts("b"));}
}
