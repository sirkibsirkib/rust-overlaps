use structs::solutions::Orientation;
use structs::solutions::*;
use structs::run_config::*;

use std::cmp::{min, max};
use std::mem::swap;
use std::collections::HashSet;

use std::thread;
use std::io::Write;
use std::io::stdout;

use bio::alignment::distance::*;

/*
only needed if reversals are enabled.
For each input string, two unique text entries are appended, one forwards and one backwards.
The orientation of these strings corresponds with the parity of their ID's
(an internal representation of the strings) with NORMAL strings having EVEN IDs and REVERSED strings having ODD parity.
Note that this is "normal" and "reversed" considering the inverse reading direction of the index.
ie: given a string "XYZ", this will append the "normal" (ZYX) to the text with ID 0 and (XYZ) to the text with ID 1.

This function simply finds the ID of the strings with the given ID such that the two
strings are from the same input string, but lie in opposite directions
*/
pub fn companion_id(id : usize) -> usize{
    if id%2==0 {id+1} else {id-1}
}

/*
Another major step in the program, the candidate verification step. (AKA the filter step)
This function returns a set of solutions, each of which corresponds to a candidate in the input set.
Only candidates that are found (somewhat naively) to have small enough error distances (as defined in config)
correspond with an output solution. Other candidates are "filtered" out.
*/
pub fn verify_all(id_a : usize, candidates : HashSet<Candidate>, config : &Config, maps : &Maps) -> HashSet<Solution>{
    let num_cands = candidates.len();
    let mut solution_set : HashSet<Solution> = HashSet::new();
    if num_cands == 0{
        return solution_set;
    }
    for c in candidates {
        let format = format!("{:#?}", &c);
        if let Some(solution) = verify(id_a, c, config, maps){
//            println!("SOL {:#?}",&solution);
            println!("{:?}\n{:#?}\n\n", format, &solution);
            solution_set.insert(solution);
        }
    }
    if config.verbose {println!("OK finished solutions for  '{}'. Verified {}/{}\t({:.2}%).",
                                maps.get_name_for(id_a), solution_set.len(), num_cands,
                                (solution_set.len() as f32) / (num_cands as f32) * 100.0)};
    solution_set
}

/*
Returns a solution corresponding with the given candidate if appropriate.
This function performs the CHECK if the candidate verifies.

The index can generate candidates that come in two forms:
>Suff-pref overlaps
a: [a1|a2]     a3==0
b:    [b2|b3]  b1==0

>Inclusions
a:    [a2]     a1==a3==0
b: [b1|b2|b3]

where a1,a2...b3 correspond with the LENGTHS of chunks of the pattern and match strings respectively,
a2 and b2 are the overlapping sections, and a1,a3,b1,b3 are the lengths of parts before and after.
*/
fn verify(id_a : usize, c : Candidate, config : &Config, maps : &Maps) -> Option<Solution>{
    let a_len = maps.get_length(id_a);
    let b_len = maps.get_length(c.id_b);
    assert_eq!(c.a3(a_len), 0);
    assert!(c.b3(b_len) >= 0); //would actually panic before this line because b3() -> usize but OK

    let a_part : &[u8] = &maps.get_string(id_a)  [c.a1()..(c.a1()+c.a2())];
    let b_part : &[u8] = &maps.get_string(c.id_b)[c.b1()..(c.b1()+c.b2())];

    println!("{:?}\n{:?}\n\n", String::from_utf8_lossy(a_part), String::from_utf8_lossy(b_part));

    let errors : u32 = if config.edit_distance{
        levenshtein(a_part, b_part)
    }else{
        assert!(a_part.len() == b_part.len());
        hamming(a_part, b_part) as u32
    };
    let mut cigar = String::new();
    cigar.push_str(&format!("ids: {}->{}", id_a, c.id_b));
    let k_limit = (config.err_rate*(max(c.overlap_a, c.overlap_b) as f32)) as u32;
    if errors <= k_limit{
        Some(solution_from_candidate(c, id_a, cigar, errors, maps, config))
    }else{
        None
    }
}

/*
Translates the input Candidate to a Solution.
This function does NOT check whether the input candidate is for a real solution.

This is one of the most confusing parts of the entire program, as everything is reversed
several times and it gets hard to keep track of how many times something is flipped.
Solutions correspond exactly with the EXTERNAL representations of the input strings,
but Candidates are largely INTERNAL (as verifying them requires the use of the index text).

*See annotation for verify() above for an explanation of a1,a2,a3,b1,b2,b3 etc. used here.
*/
fn solution_from_candidate(c : Candidate, mut id_a : usize, cigar : String, errors : u32,
                           maps : &Maps, config : &Config) -> Solution {
    let mut a_len = maps.get_length(id_a);
    let mut b_len = maps.get_length(c.id_b);
    let orientation = if !config.reversals || c.id_b%2==0{
        Orientation::Normal
    }else{
        Orientation::Reversed
    };
    let mut sol = Solution{
        id_a : id_a,
        id_b : c.id_b,
        orientation : orientation,
        overlap_a : c.overlap_a,
        overlap_b : c.overlap_b,
        overhang_left_a : c.overhang_left_a,
        overhang_right_b : (c.b3(b_len) as i32) - (c.a3(a_len) as i32),
        errors : errors,
        cigar : cigar,
    };
    if id_a > c.id_b {
        //flip which strings are A and B, as the output needs to adhere to an ascending ordering
        sol.v_flip();
    }
    sol.un_reverse(); //finally, compensate for the index being entirely backwards
    assert!(!config.reversals || sol.id_a % 2 == 0);
    sol
}