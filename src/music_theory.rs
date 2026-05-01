//! Music theory: frequencies, intervals, scales, tuning.

use crate::powers;use crate::exponential;use crate::constants::PI;

/// MIDI note to frequency (A4=440Hz).
pub fn midi_to_freq(note:u8)->f64{440.0*powers::pow(2.0,(note as f64-69.0)/12.0)}
/// Frequency to MIDI note.
pub fn freq_to_midi(freq:f64)->u8{(69.0+12.0*exponential::log2(freq/440.0)).round()as u8}

/// Note name from MIDI number.
pub fn midi_to_name(note:u8)->String{
    let names=["C","C#","D","D#","E","F","F#","G","G#","A","A#","B"];
    let octave=(note as i32/12)-1;format!("{}{}",names[(note%12)as usize],octave)
}

/// Interval ratio in equal temperament.
pub fn equal_temperament_ratio(semitones:f64)->f64{powers::pow(2.0,semitones/12.0)}

/// Just intonation ratios for common intervals.
pub fn just_ratio(interval:&str)->f64{
    match interval{"unison"=>1.0,"minor_second"=>16.0/15.0,"major_second"=>9.0/8.0,
        "minor_third"=>6.0/5.0,"major_third"=>5.0/4.0,"perfect_fourth"=>4.0/3.0,
        "tritone"=>45.0/32.0,"perfect_fifth"=>3.0/2.0,"minor_sixth"=>8.0/5.0,
        "major_sixth"=>5.0/3.0,"minor_seventh"=>9.0/5.0,"major_seventh"=>15.0/8.0,
        "octave"=>2.0,_=>1.0}
}

/// Generate major scale frequencies from root frequency.
pub fn major_scale(root:f64)->Vec<f64>{
    let steps=[0,2,4,5,7,9,11,12];
    steps.iter().map(|&s|root*equal_temperament_ratio(s as f64)).collect()
}

/// Minor (natural) scale.
pub fn minor_scale(root:f64)->Vec<f64>{
    let steps=[0,2,3,5,7,8,10,12];
    steps.iter().map(|&s|root*equal_temperament_ratio(s as f64)).collect()
}

/// Pentatonic scale.
pub fn pentatonic_scale(root:f64)->Vec<f64>{
    let steps=[0,2,4,7,9,12];
    steps.iter().map(|&s|root*equal_temperament_ratio(s as f64)).collect()
}

/// Blues scale.
pub fn blues_scale(root:f64)->Vec<f64>{
    let steps=[0,3,5,6,7,10,12];
    steps.iter().map(|&s|root*equal_temperament_ratio(s as f64)).collect()
}

/// Chord frequencies: major triad.
pub fn major_chord(root:f64)->Vec<f64>{vec![root,root*equal_temperament_ratio(4.0),root*equal_temperament_ratio(7.0)]}
/// Minor triad.
pub fn minor_chord(root:f64)->Vec<f64>{vec![root,root*equal_temperament_ratio(3.0),root*equal_temperament_ratio(7.0)]}
/// Dominant 7th.
pub fn dom7_chord(root:f64)->Vec<f64>{vec![root,root*equal_temperament_ratio(4.0),root*equal_temperament_ratio(7.0),root*equal_temperament_ratio(10.0)]}

/// Cents between two frequencies.
pub fn cents(f1:f64,f2:f64)->f64{1200.0*exponential::log2(f2/f1)}

/// Frequency of nth harmonic.
pub fn harmonic(fundamental:f64,n:u32)->f64{fundamental*n as f64}

/// Beat frequency between two tones.
pub fn beat_frequency(f1:f64,f2:f64)->f64{crate::arithmetic::abs(f1-f2)}

/// BPM to period (seconds per beat).
pub fn bpm_to_period(bpm:f64)->f64{60.0/bpm}

/// Note duration in seconds (given BPM and beat fraction).
pub fn note_duration(bpm:f64,fraction:f64)->f64{60.0/bpm*fraction}

/// Pythagorean tuning: generate scale by stacking perfect fifths (3/2).
pub fn pythagorean_scale(root:f64)->Vec<f64>{
    let mut notes:Vec<f64>=(0..12).map(|i|{let mut f=root;
        for _ in 0..i{f*=3.0/2.0;while f>=root*2.0{f/=2.0;}}f}).collect();
    notes.sort_by(|a,b|a.partial_cmp(b).unwrap());notes
}

/// Consonance score (lower = more consonant).
pub fn consonance(f1:f64,f2:f64)->f64{
    let ratio=if f1<f2{f2/f1}else{f1/f2};
    // Approximate: check closeness to simple ratios
    let simple=[(1.0,1),(2.0,1),(3.0,2),(4.0,3),(5.0,4),(5.0,3),(6.0,5)];
    simple.iter().map(|&(n,d)|crate::arithmetic::abs(ratio-n/d as f64)).fold(f64::INFINITY,f64::min)
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_a4(){assert!(approx_eq(midi_to_freq(69),440.0,1e-6));}
    #[test]fn test_freq_to_midi(){assert_eq!(freq_to_midi(440.0),69);}
    #[test]fn test_midi_name(){assert_eq!(midi_to_name(60),"C4");}
    #[test]fn test_octave(){assert!(approx_eq(equal_temperament_ratio(12.0),2.0,1e-10));}
    #[test]fn test_cents(){assert!(approx_eq(cents(440.0,880.0),1200.0,1e-6));}
    #[test]fn test_major(){let s=major_scale(440.0);assert_eq!(s.len(),8);}
}
