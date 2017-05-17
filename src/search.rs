
use bio::data_structures::fmindex::*;
use bio::data_structures::suffix_array::RawSuffixArray;
use bio::data_structures::fmindex::FMIndexable;

use std::collections::HashSet;
use std;
use std::cmp::max;


use structs::run_config::*;
use structs::solutions::*;

use algorithm_modes::kucherov::get_block_lengths;
use algorithm_modes::kucherov::candidate_condition;
use algorithm_modes::kucherov::filter_func;


pub static ALPH : &'static [u8] = b"ACGNT";
pub static READ_ERR : u8 = b'N';

pub trait GeneratesCandidates : FMIndexable {

    fn generate_candidates(&self,
                           pattern : &[u8],
                           config : &Config,
                           maps : &Maps,
                           id_a : usize,
                           sa : &RawSuffixArray,
                ) -> HashSet<Candidate> {

        println!("WORKING FOR id_a == {}", id_a);
        println!("\tPATTERN {}", String::from_utf8_lossy(pattern));
        let patt_len = pattern.len();
        let block_lengths = get_block_lengths(patt_len as i32, config.err_rate, config.thresh);
        let mut candidate_set: HashSet<Candidate> = HashSet::new();
        let max_b_len =
            if config.reversals {patt_len} else {(patt_len as f32 / (1.0 - config.err_rate)).floor() as usize};
        let block_id_lookup = get_block_id_lookup(max_b_len, config, &block_lengths);
        println!("block_id_lookup {:?}", &block_id_lookup);
        let full_interval = Interval {
            lower: 0,
            upper: self.bwt().len() - 1,
        };

        let mut p_i : i32 = (patt_len-1) as i32;

        // FOR EACH FILTER, from smallest prefix to total prefix
        // PREFIX FILTERS
        println!("\n\nGenerate candidates!");
        let patt_blocks : i32 = block_lengths.len() as i32;
        for (block_id, block_len) in block_lengths.iter().enumerate() {
            println!("\n\n\nnew suff of len {} for {}", block_id+1, String::from_utf8_lossy(&pattern));
            let cns = SearchConstants{
                pattern: pattern,
                config : config,
                maps : maps,
                block_id_lookup : &block_id_lookup,
                sa : sa,
                id_a : id_a,
                first_block_id : block_id as i32,
                patt_blocks : patt_blocks,
                blind_chars : patt_len - p_i as usize - 1,
            };
            println!(">> id_a {} first_block_id {} blind chars {}",
                     cns.id_a, cns.first_block_id, cns.blind_chars);
            self.recurse_candidates(&mut candidate_set, &cns, 0, p_i, 0, 0, 0, &full_interval, &String::new());
            p_i -= *block_len;
        }
        println!("DONE FOR id_a == {}", id_a);
        candidate_set
    }


    fn recurse_candidates(&self,
                          cand_set : &mut HashSet<Candidate>,
                          cns : &SearchConstants,
                          errors : i32,
                          p_i : i32,
                          indel_balance : i32,
                          a_match_len : usize,
                          b_match_len : usize,
                          matches : &Interval,
                          debug : &str){
        println!("'{}'   recurse p_i=={}    {:?}", debug, p_i, &matches);
        if matches.lower > matches.upper{
            // range is (lower, upper)  ie. both inclusive
            return
        }

        let completed_blocks : i32 = cns.block_id_lookup.get(max(0, p_i) as usize).expect("COMPLETED BLOCKS") - cns.first_block_id;
        let permitted_error : i32 =
            filter_func(completed_blocks, cns.patt_blocks);
        //        let has_spare_error : bool = errors < permitted_error;

        // ADD SUFFIX CANDIDATES
        let generous_match_len = std::cmp::max(a_match_len, b_match_len);
        let cand_condition_satisfied =
            candidate_condition(generous_match_len as i32, completed_blocks, cns.config.thresh, errors);
        if cand_condition_satisfied {
            println!("  !!!!!! CAND COND SATISFIED! :D {}", debug);
            let a = b'$';
            let less = self.less(a);
            let dollar_interval = Interval {
                lower : less + if matches.lower > 0 { self.occ(matches.lower - 1, a) } else { 0 },
                upper : less + self.occ(matches.upper, a),
            };
            //            println!("interval {:?}    ==$==> {:?}", &matches, &dollar_interval);
            let positions = dollar_interval.occ(cns.sa);
            for p in positions {
                let fetch_index = (p + 1) as usize;
                println!("FETCHING ID FOR INDEX {}", fetch_index);

                let id_b = *cns.maps.id2index_bdmap.get_by_second(&fetch_index).expect("DOLLAR MAP BAD");
                let c = Candidate {
                    id_b: id_b,
                    overlap_a: a_match_len + cns.blind_chars,
                    overlap_b: b_match_len + cns.blind_chars,
                    overhang_right_b: 0,
                    debug_str : debug.to_owned(),
                };
                if id_b == cns.id_a{
                    //skip obvious self-matches
                    println!("  !!!!!! SELF CAND {:?}, {:?}", p, debug);
                    cand_set.insert(c);
                    continue;
                }

                println!("  !!!!!! ~~~  adding FOUND CANDIDATE AT {} with {}", p, debug);
                cand_set.insert(c);
            }
        }
        if p_i == -1{
            //END OF PATTERN
            //INCLUSIONS GO HERE
            println!("KILL");
            return
        }

        // walking + substitution match
        for &a in ALPH.iter(){
            let less = self.less(a);
            //            println!("less  {}", less);
            //            println!("match is {:?}", matches);
            let next_interval = Interval{
                lower : less + if matches.lower > 0 { self.occ(matches.lower - 1, a) } else { 0 },
                upper : less + self.occ(matches.upper, a) - 1,
            };
            let p_char = *cns.pattern.get(p_i as usize).expect("THE P CHAR");
            let recurse_error =  if p_char == a && a != READ_ERR {errors} else {errors + 1};
            if recurse_error <= permitted_error{
                let mut next_debug = format!("{}{}", a as char, debug);
                self.recurse_candidates(cand_set,
                                        cns,
                                        recurse_error,
                                        p_i-1,
                                        0,
                                        a_match_len+1,
                                        b_match_len+1,
                                        &next_interval,
                                        &next_debug);
            }
        }
    }
}


#[derive(Debug)]
struct SearchConstants<'a>{
    config : &'a Config,
    maps : &'a Maps,

    block_id_lookup : &'a Vec<i32>,
    sa : &'a RawSuffixArray,
    pattern: &'a [u8],
    id_a : usize,
    blind_chars : usize,

    first_block_id : i32,
    patt_blocks : i32,
}

fn get_block_id_lookup(max_b_len : usize, config : &Config, block_lengths : &[i32]) -> Vec<i32>{
    let mut lookup : Vec<i32> = Vec::new();
    for (id, block_length) in block_lengths.iter().enumerate() {
        for _ in 0..*block_length{
            lookup.push(id as i32);
        }
    }
    lookup.reverse();
    lookup.shrink_to_fit();
    lookup
}