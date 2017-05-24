use std::cmp::max;
use std::collections::HashSet;

use bio::alignment::distance::*;

use structs::solutions::*;
use structs::run_config::*;

use useful::*;


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
//        let format = format!("{:#?}", &c);
        if let Some(solution) = verify(id_a, c, config, maps){
//            println!("SOL {:#?}",&solution);
//            println!("{:?}\n{:#?}\n\n", format, &solution);
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
    assert_eq!(c.a3(a_len), 0);
    //b3 is usize, so implicitly b3 >= 0

    let a_part : &[u8] = &maps.get_string(id_a)  [c.a1()..(c.a1()+c.a2())];
    let b_part : &[u8] = &maps.get_string(c.id_b)[c.b1()..(c.b1()+c.b2())];

    let errors : u32 = if config.edit_distance{
        levenshtein(a_part, b_part)
    }else{
        assert!(a_part.len() == b_part.len());
        hamming(a_part, b_part) as u32
    };
//    let mut cigar = String::new();
//    cigar.push_str(&format!("ids: {}->{}", id_a, c.id_b));
    let k_limit = (config.err_rate*(max(c.overlap_a, c.overlap_b) as f32)) as u32;
    if errors <= k_limit{
        Some(solution_from_candidate(c, id_a, errors, maps, config))
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
fn solution_from_candidate(c : Candidate, id_a : usize, errors : u32,
                           maps : &Maps, config : &Config) -> Solution {
    let a_len = maps.get_length(id_a);
    let b_len = maps.get_length(c.id_b);
    let orientation = relative_orientation(id_a, c.id_b, config.reversals);
    let mut sol = Solution{
        id_a : id_a,
        id_b : c.id_b,
        orientation : orientation,
        overlap_a : c.overlap_a,
        overlap_b : c.overlap_b,
        overhang_left_a : c.overhang_left_a,
        overhang_right_b : (c.b3(b_len) as i32) - (c.a3(a_len) as i32),
        errors : errors,
//        cigar : cigar,
    };
    translate_solution_to_external(&mut sol, config);
    sol
}

fn translate_solution_to_external(sol : &mut Solution, config : &Config){
    assert!(sol.id_a != sol.id_b);
    if config.reversals {
        assert!(sol.id_a != companion_id(sol.id_b, config.reversals));
    }

    if sol.id_a > sol.id_b {
        sol.v_flip();
    }
    assert!(sol.id_a <= sol.id_b);

    //TODO reversals allow for things to be skipped.

    if config.reversals {
        if for_reversed_string(sol.id_a){
            sol.h_flip(config.reversals);
        }
        assert!(!for_reversed_string(sol.id_a));
    }

    sol.un_reverse(); //finally, compensate for the index being entirely backwards
    assert!(!config.reversals || sol.id_a % 2 == 0);
}