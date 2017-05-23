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

pub fn companion_id(id : usize) -> usize{
    if id%2==0 {id+1} else {id-1}
}

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

fn verify(id_a : usize, c : Candidate, config : &Config, maps : &Maps) -> Option<Solution>{
    let a_len = maps.get_length(id_a);
    let b_len = maps.get_length(c.id_b);
    assert_eq!(c.a3(a_len), 0);
    assert!(c.b3(b_len) >= 0); //would actually crash because usize but OK

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
        sol.v_flip();
    }
    sol.un_reverse();
    assert!(!config.reversals || sol.id_a % 2 == 0);
    sol


//    let mut id_b = c.id_b;

//    let mut overlap_a = c.overlap_a;
//    let mut overlap_b = c.overlap_b;
//    let mut overhang_left_a = c.overhang_left_a;
//    let mut overhang_right_b = if overhang_left_a == 0{
//        (b_len as i32) - (overlap_b as i32)
////        (b_len as i32) - (a_len as i32) - overhang_left_a;
//    } else {
//        (b_len as i32) - (overlap_b as i32) - overhang_left_a
//    };
//    println!("!??!?  {}", overhang_right_b);


    // GUARANTEE 1/2: id_a <= id_b
    // REMEDY: vertical flip. a becomes b, b becomes a.
//    if id_a > id_b {
//        //
////        println!("VFLIP");
//        overhang_left_a *= -1;
//        overhang_right_b *= -1;
//        swap(&mut a_len, &mut b_len);
//        swap(&mut id_a, &mut id_b);
//        swap(&mut overlap_a, &mut overlap_b);
////        cigar = cigar.vflip();  // I->D, D->I
//    }
//
////    // GUARANTEE 2/2: A is a string from the input (not a flipped string)
////    // REMEDY: horizontal flip. a becomes b, b becomes a.
////    if config.reversals && id_a % 2 == 1 {
//////        println!("HFLIP");
////        swap(&mut overhang_left_a, &mut overhang_right_b);
////        overhang_left_a *= -1;
////        overhang_right_b *= -1;
////        id_a = companion_id(id_a);
////        id_b = companion_id(id_b);
//////        cigar.h_flip(); // XYZ -> ZYX
////    }
//    let orientation = if !config.reversals || id_b%2==0{
//        Orientation::Normal
//    }else{
//        Orientation::Reversed
//    };
//
//    // INTERNAL --> EXTERNAL REPRESENTATION
//    swap(&mut overhang_left_a, &mut overhang_right_b);
//    overhang_left_a *= -1;
//    overhang_right_b *= -1;
//
//
//    Solution{
//        id_a : id_a,
//        id_b : id_b,
//        orientation : orientation,
//        overlap_a : overlap_a,
//        overlap_b : overlap_b,
//        overhang_left_a : overhang_left_a,
//        overhang_right_b : overhang_right_b,
//        errors : errors,
//        cigar : cigar,
//    }
}