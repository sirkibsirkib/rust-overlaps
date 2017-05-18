use structs::solutions::Orientation;
use structs::solutions::*;
use structs::run_config::*;

use std::cmp::max;
use std::mem::swap;
use std::collections::HashSet;

use bio::alignment::distance::*;

pub fn companion_id(id : usize) -> usize{
    if id%2==0 {id+1} else {id-1}
}

pub fn verify_all(id_a : usize, candidates : HashSet<Candidate>, config : &Config, maps : &Maps) -> HashSet<Solution>{
    let mut solution_set : HashSet<Solution> = HashSet::new();
    for c in candidates {
        if let Some(solution) = verify(id_a, c, config, maps){
            solution_set.insert(solution);
        }
    }
    if config.verbose {println!("OK finished solutions for  '{}'.", maps.get_name_for(id_a))};
    solution_set
}

fn verify(id_a : usize, c : Candidate, config : &Config, maps : &Maps) -> Option<Solution>{

    let a_len = maps.get_length(id_a);
    let b_len = maps.get_length(c.id_b);

    // the search guarantees this; Thus omitted from the candidate struct.
    let overhang_left_a : i32 = a_len as i32 - c.overlap_a as i32 ;
    let overlap_a_start : usize = max(0, overhang_left_a) as usize;
    let overlap_a_end : usize = overlap_a_start + c.overlap_a;

    if b_len < c.overlap_b{
        // the blind spot goes over the end of B string. (bad match!)
        return None;
    }
    let overlap_b_start : usize  = max(0, -overhang_left_a) as usize;
    let overlap_b_end : usize  = overlap_b_start + c.overlap_b as usize;

    let a_part : &[u8] = &maps.get_string(id_a)  [overlap_a_start..overlap_a_end];
    let b_part : &[u8] = &maps.get_string(c.id_b)[overlap_b_start..overlap_b_end];

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
        Some(solution_from_candidate(c, id_a, cigar, errors, maps, config, overhang_left_a))
    }else{
        None
    }
}

fn solution_from_candidate(c : Candidate, mut id_a : usize, cigar : String, errors : u32,
                           maps : &Maps, config : &Config, mut overhang_left_a : i32) -> Solution {
    let mut a_len = maps.get_length(id_a);
    let mut b_len = maps.get_length(c.id_b);
    let mut id_b = c.id_b;

    let mut overlap_a = c.overlap_a;
    let mut overlap_b = c.overlap_b;
    let mut overhang_right_b = c.overhang_right_b;

    assert!(overhang_right_b != 8 && overhang_right_b != -8);
    assert!(overhang_left_a != 8 && overhang_left_a != -8);

    // GUARANTEE 1/2: id_a <= id_b
    // REMEDY: vertical flip. a becomes b, b becomes a.
    if id_a > id_b {
        //
        overhang_left_a *= -1;
        overhang_right_b *= -1;
        swap(&mut a_len, &mut b_len);
        swap(&mut id_a, &mut id_b);
        swap(&mut overlap_a, &mut overlap_b);
//        cigar = cigar.vflip();  // I->D, D->I
    }

    // GUARANTEE 2/2: A is a string from the input (not a flipped string)
    // REMEDY: horizontal flip. a becomes b, b becomes a.
    if config.reversals && id_a % 2 == 1 {
        swap(&mut overhang_left_a, &mut overhang_right_b);
        overhang_left_a *= -1;
        overhang_right_b *= -1;
        id_a = companion_id(id_a);
        id_b = companion_id(id_b);
//        cigar.h_flip(); // XYZ -> ZYX
    }
    let orientation = if !config.reversals || id_b%2==0{
        Orientation::Normal
    }else{
        Orientation::Reversed
    };

    Solution{
        id_a : id_a,
        id_b : id_b,
        orientation : orientation,
        overlap_a : overlap_a,
        overlap_b : overlap_b,
        overhang_left_a : overhang_left_a,
        overhang_right_b : overhang_right_b,
        errors : errors,
        cigar : cigar,
    }
}