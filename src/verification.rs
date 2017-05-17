use std::{thread, time};
use std::collections::HashSet;

use structs::solutions::Orientation;
use structs::solutions::*;
use structs::run_config::*;

use std::cmp::{max, min};
use std::mem::swap;

use bio::alignment::distance::*;

fn sleep(){
    thread::sleep(time::Duration::from_millis(1000));
}

fn companion_id(id : usize) -> usize{
    if id%2==0 {id+1} else {id-1}
}

pub fn verify_all(id_a : usize, candidates : HashSet<Candidate>, config : &Config, maps : &Maps) -> HashSet<Solution>{
    println!("VERIFY FOR ID {} START", id_a);
    let mut solution_set : HashSet<Solution> = HashSet::new();
    for c in candidates{
        println!("Cand looks like {:#?}", &c);
        match verify(id_a, c, config, maps){
            Some(x) => {solution_set.insert(x);},
            None => (),
        }
//        if let Some(x) = verify(id_a, c, config, maps){
//            println!("SUCCEEDED");
//            solution_set.insert(x);
//        }else{
//            println!("failed!");
//        }
    }
    println!("VERIFY FOR ID {} END. found {} solutions", id_a, solution_set.len());
    solution_set
}

// incoming cands are already oriented
fn verify(id_a : usize, c : Candidate, config : &Config, maps : &Maps) -> Option<Solution>{
    println!("VERIFY a/b ({}, {})", id_a, c.id_b);

    let a_len = maps.get_length(id_a);
    let b_len = maps.get_length(c.id_b);
    let overhang_left_a : i32 = (a_len - c.overlap_a) as i32;

    let a_start : usize = max(0, overhang_left_a) as usize;
    let a_end : usize = a_start + c.overlap_a;

    //[pattern str]
    //       [match str]
    // but might be
    //[pattern~~~~~]
    //     [match?????  where pattern overlap extends over the end
    // check if

    println!("===============");
    for _ in 0..-overhang_left_a{print!{" "};}
    println!("{}", String::from_utf8_lossy(maps.get_string(id_a)));
    for _ in 0..overhang_left_a{print!{" "};}
    println!("{}", String::from_utf8_lossy(maps.get_string(c.id_b)));
    println!("===============");

    return None;
    if b_len < c.overlap_b{
        println!("HANGING OVER THE EDGE");
        return None;
    }


//    let a_end2 : usize  = min(a_len, a_len-(c.overhang_right_b as usize));
//    //TODO return None here if you detect you're out of bounds of one of the strings
//    assert_eq!(a_end, a_end2);          //DEBUG! can remove a_end2 if this checks out

//    println!("SO FAR SO GOOD");
//    return None;
    let b_start : usize  = max(0, -overhang_left_a) as usize;
    let b_end : usize  = b_start + c.overlap_b as usize;
//    let b_end2 : usize  = min(b_len, b_len+c.overhang_right_b as usize);
//    //TODO return None here if you detect you're out of bounds of one of the strings
//    assert_eq!(b_end, b_end2);          //DEBUG! can remove b_end2 if this checks out


    let a_part : &[u8] = &maps.get_string(id_a)  [a_start..a_end];
    let b_part : &[u8] = &maps.get_string(c.id_b)[b_start..b_end];

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

fn solution_from_candidate(c : Candidate, mut id_a : usize, mut cigar : String, errors : u32, maps : &Maps, config : &Config) -> Solution {
    let mut a_len = maps.get_length(id_a);
    let mut b_len = maps.get_length(c.id_b);
    let mut id_b = c.id_b;

    let mut overlap_a = c.overlap_a;
    let mut overlap_b = c.overlap_b;
    let (mut overhang_left_a, mut overhang_right_b) = if c.overhang_right_b == 0 {
        (
            -((a_len - overlap_a) as i32),
            -((b_len - overlap_b) as i32),
        )
    } else {
        (
            -(a_len as i32 + c.overhang_right_b - overlap_a as i32),
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


    println!("  ==>");


    for _ in 0..-overhang_left_a{print!{" "};}
    println!("{}", String::from_utf8_lossy(maps.get_string(id_a)));
    for _ in 0..overhang_left_a{print!{" "};}
    println!("{}", String::from_utf8_lossy(maps.get_string(c.id_b)));

    let orientation = if !config.reversals || id_b%2==0{
        Orientation::Normal
    }else{
        Orientation::Reversed
    };
    let a_name = maps.id2name_vec.get(id_a);
    let b_name = maps.id2name_vec.get(id_b);

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