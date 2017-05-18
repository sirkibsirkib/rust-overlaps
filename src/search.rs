
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
use std::io::stdout;
use std::io::Write;


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

//        println!("WORKING FOR id_a == {}", id_a);
//        println!("\tPATTERN {}", String::from_utf8_lossy(pattern));
        let patt_len = pattern.len();
        let block_lengths = get_block_lengths(patt_len as i32, config.err_rate, config.thresh);
        let mut candidate_set: HashSet<Candidate> = HashSet::new();
        let max_b_len =
            if config.reversals {patt_len} else {(patt_len as f32 / (1.0 - config.err_rate)).floor() as usize};
        let block_id_lookup = get_block_id_lookup(max_b_len, config, &block_lengths);
//        println!("block_id_lookup {:?}", &block_id_lookup);
        let full_interval = Interval {
            lower: 0,
            upper: self.bwt().len() - 1,
        };

        let mut p_i : i32 = (patt_len-1) as i32;

        // FOR EACH FILTER, from smallest prefix to total prefix
        // PREFIX FILTERS
//        println!("\n\nGenerate candidates! for ID {}", id_a);
        let patt_blocks : i32 = block_lengths.len() as i32;
        for (block_id, block_len) in block_lengths.iter().enumerate() {
//            println!("\nnew suff of len {} for {}", block_id+1, String::from_utf8_lossy(&pattern));
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
//            println!(">> id_a {} first_block_id {} blind chars {}",
//                     cns.id_a, cns.first_block_id, cns.blind_chars);
            self.recurse_candidates(&mut candidate_set, &cns, 0, p_i, 0, 0, 0, &full_interval, &String::new());
            p_i -= *block_len;
        }
//        println!("DONE FOR id_a == {}. found {} unique cands", id_a, candidate_set.len());
        if config.verbose {println!("OK finished candidates for '{}'.", maps.id2name_vec.get(id_a).expect("UNKNOWN ID"))};
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
                          match_interval : &Interval,
                          debug : &str){
        if match_interval.lower > match_interval.upper{
            // range is (lower, upper)  ie. both inclusive
            return
        }

        let completed_blocks : i32 = match cns.block_id_lookup.get(p_i as usize){
            Some(x) => x - cns.first_block_id,
            None    => cns.patt_blocks - cns.first_block_id,
        };

        let permitted_errors : i32 =
            filter_func(completed_blocks, cns.patt_blocks);
        //        let has_spare_error : bool = errors < permitted_error;

        // ADD SUFFIX CANDIDATES
        let generous_match_len = std::cmp::max(a_match_len, b_match_len) + 1;
        let cand_condition_satisfied =
            candidate_condition(generous_match_len as i32, completed_blocks, cns.config.thresh, errors);



//        println!("'{}'{}   recurse p_i=={}    {:?}", debug, if cand_condition_satisfied {'*'} else {' '},p_i, &match_interval);
        if cand_condition_satisfied {
//            println!("CAND? YES");
            let a = b'$';
            let less = self.less(a);
            let dollar_interval = Interval {
                lower : less + if match_interval.lower > 0 { self.occ(match_interval.lower - 1, a) } else { 0 },
                upper : less + self.occ(match_interval.upper, a),
            };
            //            println!("interval {:?}    ==$==> {:?}", &matches, &dollar_interval);
            let positions = dollar_interval.occ(cns.sa);
            add_candidate_here(positions, cand_set, cns, a_match_len, b_match_len, debug, false);
        }

//        // insertion
//        for &a in ALPH.iter(){
//            let less = self.less(a);
//            //            println!("less  {}", less);
//            //            println!("match is {:?}", matches);
//            let next_interval = Interval{
//                lower : less + if match_interval.lower > 0 { self.occ(match_interval.lower - 1, a) } else { 0 },
//                upper : less + self.occ(match_interval.upper, a) - 1,
//            };
//            let p_char = *cns.pattern.get(p_i as usize).expect("THE P CHAR");
//            let recurse_error =  if p_char == a && a != READ_ERR {errors} else {errors + 1};
//            if recurse_error <= permitted_error{
//                let mut next_debug = format!("{}{}", a as char, debug);
//                self.recurse_candidates(cand_set,
//                                        cns,
//                                        recurse_error,
//                                        p_i-1,  //step left
//                                        0,      //indel balance reset
//                                        a_match_len + 1,
//                                        b_match_len + 1,
//                                        &next_interval,
//                                        &next_debug);
//            }
//        }

        if p_i == -1{
            if cns.config.inclusions && cand_condition_satisfied{
                let inclusion_interval = Interval{
                    lower : match_interval.lower,
                    upper : match_interval.upper + 1,
                }; // final interval must have exclusive end
                let positions = inclusion_interval.occ(cns.sa);
//                println!("{} cands?", positions.len());
                add_candidate_here(positions, cand_set, cns, a_match_len, b_match_len, debug, true);
            }
            if cns.config.edit_distance{
                //Nothing do do here if we are not doing edit distance
                return;
            }
        }

        for &a in ALPH.iter(){
            let less = self.less(a);
            //            println!("less  {}", less);
            //            println!("match is {:?}", matches);
            let next_interval = Interval{
                lower : less + if match_interval.lower > 0 { self.occ(match_interval.lower - 1, a) } else { 0 },
                upper : less + self.occ(match_interval.upper, a) - 1,
            };
            let p_char = *cns.pattern.get(p_i as usize).expect("THE P CHAR");
            let recurse_errors =  if p_char == a && a != READ_ERR {errors} else {errors + 1};
            let mut next_debug = format!("{}{}", a as char, debug);

            //substitutions and matches
            if recurse_errors <= permitted_errors && p_i > -1{
                self.recurse_candidates(cand_set,
                                        cns,
                                        recurse_errors,
                                        p_i-1,  //step left
                                        0,      //indel balance reset
                                        a_match_len + 1,
                                        b_match_len + 1,
                                        &next_interval,
                                        &next_debug);
            }
            //insertions
            if errors < permitted_errors && cns.config.edit_distance {
                self.recurse_candidates(cand_set,
                                        cns,
                                        errors + 1, //always induces an error
                                        p_i-1,      //step left
                                        1,          //insertion balance
                                        a_match_len,//the pattern string doesn't grow
                                        b_match_len + 1,
                                        &next_interval,
                                        &next_debug);
            }
        }

        // deletion
        if cns.config.edit_distance && errors < permitted_errors{
            if indel_balance <= 0 && p_i > 1{
                let mut next_debug = format!("{}{}", '_', debug);
                self.recurse_candidates(cand_set,
                                        cns,
                                        errors + 1,
                                        p_i - 1,         //one step without matching
                                        -1,              //balance
                                        a_match_len + 1,
                                        b_match_len,     //the matched string doesn't grow
                                        &match_interval, //stays unchanged
                                        &next_debug);
            }
        }
    }
}

use verification;

#[inline]
fn add_candidate_here(positions : Vec<usize>,
                      cand_set : &mut HashSet<Candidate>,
                      cns : &SearchConstants, a_match_len : usize,
                      b_match_len : usize, debug : &str, inclusion : bool){
    for p in positions {
        let id_b = if inclusion {
            cns.maps.find_id_for_index_within(p)
        } else {
            *cns.maps.id2index_bdmap.get_by_second(&(p+1)).expect("UH OH")
        } as usize;
        if id_b == cns.id_a ||
            (cns.config.edit_distance && cns.id_a == verification::companion_id(cns.id_a)){
            //self-match or partner match
            continue;
        }
        let overlap_a = a_match_len + cns.blind_chars;
        let overlap_b = b_match_len + cns.blind_chars;
        let a_len = cns.pattern.len();
        let b_len = cns.maps.get_length(id_b);
        if overlap_a == a_len && overlap_b == b_len && cns.id_a > id_b {
            //perfect complete overlap. this one is deemed to be redundant
            continue;
        }
        let overhang_right_b = b_len as i32 - (overlap_b as i32);
        if overhang_right_b < 0{
            //Out of right-side bounds
            continue;
        }
        let c = Candidate {
            id_b: id_b,
            overlap_a: overlap_a,
            overlap_b: overlap_b,
            overhang_right_b: overhang_right_b,
            debug_str : debug.to_owned(),
        };
        cand_set.insert(c);
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