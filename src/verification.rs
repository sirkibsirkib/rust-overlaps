
use std::io;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::fs;
use std::fs::File;
use std::io::{Write, BufWriter};

use structs::solutions::Orientation;
use structs::solutions::*;
use structs::run_config::*;

use std::cmp::{max, min};
use std::mem::swap;

use bio::alignment::distance::*;


fn companion_id(id : i32) -> i32{
    if id%2==0 {id+1} else {id-1}
}

pub fn verify_all(id_a : i32, candidates : HashSet<Candidate>, config : &Config, maps : &Maps) -> HashSet<Solution>{
    let mut solution_set : HashSet<Solution> = HashSet::new();
    for c in candidates{
        if let Some(x) = verify(id_a, c, config, maps){
            solution_set.insert(x);
        }
    }
    solution_set
}

// incoming cands are already oriented
fn verify(id_a : i32, c : Candidate, config : &Config, maps : &Maps) -> Option<Solution>{
    let overhang_left_a = 0;

    let a_len = maps.id2str_in_s.get(&id_a).expect("WTF").len() as i32;
    let b_len = maps.id2str_in_s.get(&c.id_b).expect("WTF").len() as i32;

    let a_start = max(0, overhang_left_a);
    let a_end = a_start + c.overlap_a;
    let a_end2 = min(a_len, a_len-c.overhang_right_b);
    //TODO return None here if you detect you're out of bounds of one of the strings
    assert_eq!(a_end, a_end2);          //DEBUG! can remove a_end2 if this checks out

    let b_start = max(0, -overhang_left_a);
    let b_end = b_start + c.overlap_b;
    let b_end2 = min(b_len, b_len+c.overhang_right_b);
    //TODO return None here if you detect you're out of bounds of one of the strings
    assert_eq!(b_end, b_end2);          //DEBUG! can remove b_end2 if this checks out

    let a_part : &[u8] = &maps.id2str_in_s.get(&id_a).expect("OMG")  [a_start as usize..a_end as usize];
    let b_part : &[u8] = &maps.id2str_in_s.get(&c.id_b).expect("OMG")[b_start as usize..b_end as usize];

    let errors : u32 = if config.edit_distance{
        assert!(a_part.len() == b_part.len());
        hamming(a_part, b_part) as u32
    }else{
        levenshtein(a_part, b_part)
    };
    let mut cigar = String::new();
    cigar.push_str("MEMES");
    let k_limit = (config.err_rate*(max(c.overlap_a, c.overlap_b) as f32)) as u32;
    if errors <= k_limit{
        Some(solution_from_candidate(c, id_a, cigar, errors, maps, config))
    }else{
        None
    }
}

fn solution_from_candidate(c : Candidate, mut id_a : i32, mut cigar : String, errors : u32, maps : &Maps, config : &Config) -> Solution {
    let mut a_len = maps.id2str_in_s.get(&id_a).expect("WTF").len() as i32;
    let mut b_len = maps.id2str_in_s.get(&c.id_b).expect("WTF").len() as i32;
    let mut id_b = c.id_b;

    let mut overlap_a = c.overlap_a as i32;
    let mut overlap_b = c.overlap_b as i32;
//    let (overhang_left_a, overhang_right_b) = if c.overhang_right_b == 0 {
//        (
//            -(b_len - overlap_b),
//            -(a_len - overlap_a),
//        )
//    } else {
//        (
//            -(b_len - c.overhang_right_b - overlap_b),
//            c.overhang_right_b,
//        )
//    };
    let (mut overhang_left_a, mut overhang_right_b) = if c.overhang_right_b == 0 {
        (
            -(b_len - overlap_b),
            -(a_len - overlap_a),
        )
    } else {
        (
            -(b_len - c.overhang_right_b - overlap_b),
            c.overhang_right_b,
        )
    };

    //START TRANSFORM

    //ensure smaller string comes first
    //vertical flip
    if id_a > id_b {
        //
        overhang_left_a *= -1;
        overhang_right_b *= -1;
        swap(&mut a_len, &mut b_len);
        swap(&mut id_a, &mut id_b);
        swap(&mut overlap_a, &mut overlap_b);
//        cigar = cigar.vflip();  // I->D, D->I
    }

    //first string is a reverse
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
    let a_name = maps.id2name.get(&id_a);
    let b_name = maps.id2name.get(&id_b);
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