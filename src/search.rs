
use bio::data_structures::fmindex::*;
use bio::data_structures::suffix_array::RawSuffixArray;
use bio::data_structures::fmindex::FMIndexable;

use std;
use std::collections::HashSet;
use std::cmp::max;
use std::hash::Hash;

////////////////////////////////

use structs::run_config::*;
use structs::solutions::*;

use verification;

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
        let patt_len = pattern.len();
        let block_lengths = get_block_lengths(patt_len as i32, config.err_rate, config.thresh);
        let mut candidate_set: HashSet<Candidate> = HashSet::new();
        let block_id_lookup = get_block_id_lookup(&block_lengths);
        let full_interval = Interval {
            lower: 0,
            upper: self.bwt().len() - 1,
        };
        let mut p_i : i32 = (patt_len-1) as i32;
        let patt_blocks : i32 = block_lengths.len() as i32;
        for (block_id, block_len) in block_lengths.iter().enumerate() {
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
            self.recurse_candidates(
                &mut candidate_set, &cns, 0, p_i,
                LastOperation::Initial, 0, 0,
                &full_interval, &String::new());

            p_i -= *block_len;
        }
        if candidate_set.is_empty(){
            if config.verbose {println!("OK no candidates found for '{}', skipping verification.",
                                        maps.get_name_for(id_a))};
        } else {
            if config.verbose {println!("OK finished candidates for '{}'.",
                                        maps.get_name_for(id_a))};
        }
        candidate_set
    }

    fn recurse_candidates(&self,
                          cand_set : &mut HashSet<Candidate>,
                          cns : &SearchConstants,
                          errors : i32,
                          p_i : i32,
                          last_operation : LastOperation,
                          a_match_len : usize,
                          b_match_len : usize,
                          match_interval : &Interval,
                          debug : &str){
        if match_interval.lower > match_interval.upper{
            // range is inclusive on both ends within the walk.
            // empty range so prune branch
            return
        }

        let completed_blocks : i32 = match cns.block_id_lookup.get(p_i as usize){
            Some(x) => x - cns.first_block_id,
            //the final step (p_i==-1) can still insert. we can't return yet, but we are beyond the blocks
            None    => cns.patt_blocks - cns.first_block_id,
        };
        let permitted_errors : i32 = filter_func(completed_blocks, cns.patt_blocks);


        let generous_match_len = std::cmp::max(a_match_len, b_match_len) + 1;
        let cand_condition_satisfied =
            candidate_condition(generous_match_len as i32, completed_blocks, cns.config.thresh, errors);

        if cand_condition_satisfied && last_operation.allows_candidates(){
            // Add candidates to set for matched b strings preceded by '$'
            let a = b'$';
            let less = self.less(a);
            let dollar_interval = Interval {
                lower : less + if match_interval.lower > 0 { self.occ(match_interval.lower - 1, a) } else { 0 },
                upper : less + self.occ(match_interval.upper, a),
            }; //final interval must have exclusive end
            let positions = dollar_interval.occ(cns.sa);
            add_candidate_here(positions, cand_set, cns, a_match_len, b_match_len, debug, false);
        }

        let pattern_finished = p_i <= -1;
        if pattern_finished {
            // end of the pattern string
            // Add inclusion candidates to set at this position for everything in the remaining range
            if cns.config.inclusions && cand_condition_satisfied && last_operation.allows_candidates(){
                let inclusion_interval = Interval{
                    lower : match_interval.lower,
                    upper : match_interval.upper + 1,
                }; // final interval must have exclusive end
                let positions = inclusion_interval.occ(cns.sa);
                add_candidate_here(positions, cand_set, cns, a_match_len, b_match_len, debug, true);
            }
            return;
            //nothing to do here.
        }

        // consider a new derived b string match, one char longer (in front) than existing match
        for &a in ALPH.iter(){
            let less = self.less(a);
            let next_interval = Interval{
                lower : less + if match_interval.lower > 0 { self.occ(match_interval.lower - 1, a) } else { 0 },
                upper : less + self.occ(match_interval.upper, a) - 1,
            };

            //TODO remove debug stuff



            let p_char = *cns.pattern.get(p_i as usize).expect("THE P CHAR");
            let recurse_errors =  if p_char == a && a != READ_ERR {errors} else {errors + 1};
            if recurse_errors <= permitted_errors {
//                let next_debug = format!("{}{}", a as char, debug);
                // recursively explore SUBSTITUTION cases (both hamming and levenshtein)
                self.recurse_candidates(cand_set,
                                        cns,
                                        recurse_errors,
                                        p_i-1,  //step left
                                        LastOperation::Substitution,
                                        a_match_len + 1,
                                        b_match_len + 1,
                                        &next_interval,
                                        debug);
            }
            if (p_char != a) && (errors < permitted_errors) && cns.config.edit_distance && last_operation.allows_insertion() {
                // recursively explore INSERTION cases (if levenshtein)
//                let next_debug = format!("{}{}", smaller(a), debug);
                self.recurse_candidates(cand_set,
                                        cns,
                                        errors + 1, //always induces an error
                                        p_i,        //don't step left
                                        LastOperation::Insertion,
                                        a_match_len,//the pattern string doesn't grow
                                        b_match_len + 1,
                                        &next_interval,
                                        debug);

            }
        }

        if cns.config.edit_distance && errors < permitted_errors && !pattern_finished{
            // recursively explore DELETION cases (if levenshtein) and have at least 1 spare pattern char to jump over
            if last_operation.allows_deletion(){

//                let next_debug = format!("{}{}", '_', debug);
                self.recurse_candidates(cand_set,
                                        cns,
                                        errors + 1,
                                        p_i - 1,         //one step without matching
                                        LastOperation::Deletion,
                                        a_match_len + 1,
                                        b_match_len,     //the matched string doesn't grow
                                        &match_interval, //stays unchanged
                                        debug);
            }
        }
    }
}

fn smaller(a : u8) -> char{
    match a as char {
        'A' => 'a',
        'C' => 'c',
        'N' => 'n',
        'G' => 'g',
        'T' => 't',
        _ => '?',
    }
}

#[derive(PartialEq, Copy, Clone)]
pub enum LastOperation{
    Initial,
    Substitution,
    Insertion,
    Deletion,
}

//fn branch_allowed(last_operation : LastOperation, desired_operation : LastOperation) -> bool{
//    match last_operation{
//        Initial => true
//    }
//}

impl LastOperation{
    fn allows_deletion(self) -> bool{
        self == LastOperation::Deletion
        || self == LastOperation::Substitution
    }

    fn allows_insertion(self) -> bool{
        self == LastOperation::Insertion
            || self == LastOperation::Substitution
    }

    fn allows_candidates(self) -> bool{
        self == LastOperation::Initial
            || self == LastOperation::Substitution
    }
}

#[inline]
fn add_candidate_here(positions : Vec<usize>,
                      cand_set : &mut HashSet<Candidate>,
                      cns : &SearchConstants, a_match_len : usize,
                      b_match_len : usize, debug : &str, inclusion : bool){
    for mut position in positions {
        if !inclusion{
            //non-inclusions include the preceding dollar sign
            position += 1;
        }

//        println!("\n\nincl: {}", inclusion);
//        cns.maps.print_text_debug();
//        println!("{}", cns.maps.push_string(debug, " ", position));

        let (id_b, index_b) = if inclusion {
            cns.maps.find_occurrence_containing(position)
        } else {
            (cns.maps.id_for(position), position)
        };
//        println!("{}", cns.maps.push_string("^", " ", index_b));
//        println!("ida {}\nidb {} indb {}", cns.id_a, id_b, index_b);


        if id_b == cns.id_a ||
            (cns.config.edit_distance &&
                cns.id_a == verification::companion_id(cns.id_a)){
            //self-match or partner match
//            println!("self match");
            continue;
        }
        let overlap_a = a_match_len + cns.blind_chars;
        let overlap_b = b_match_len + cns.blind_chars;
        let a_len = cns.pattern.len();
        let b_len = cns.maps.get_length(id_b);
        let overhang_left_a : i32 = if inclusion{
            //        [aaaa|aa]
            //   [bbbbbbbbb?????
            //   ^    ^
            //   <--->
            - (position as i32 - index_b as i32)
        } else {
            // [aaaaaa|aa]
            // <--->[b?????????
            a_len as i32 - overlap_a as i32
        };
//        println!("overhang_left_a {}", overhang_left_a);

        if overlap_a == a_len && overlap_b == b_len && cns.id_a > id_b {
            //perfect complete overlap. this one is deemed to be redundant
//            println!("redundant perfect overlap");
            continue;
        }

        if max(overhang_left_a, 0) as usize + overlap_b > b_len {
            //we have matched with a string beyond the bounds of b
//            println!("out of bounds!");
            continue;
        }

//        println!("LOHA {}", overhang_left_a);

        let mut new_debug = debug.to_owned();
        new_debug.push_str(&format!(" incl {} blind {}", inclusion, cns.blind_chars));
        let c = Candidate {
            id_b: id_b,
            overlap_a: overlap_a,
            overlap_b: overlap_b,
            overhang_left_a: overhang_left_a,
            debug_str : new_debug,
        };

        if cns.config.edit_distance{
            panic!("Problem: the amount of 'blind' characters that should be matched is essentially unknowable \
            so how t.f. am I supposed to reconcile that?");
        }
        if cand_set.contains(&c){
            println!("OLD {:#?}\nNEW {:#?}", c, c);
        }
        cand_set.insert(c);
    }
}


#[derive(Debug)]
pub struct SearchConstants<'a>{
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

fn get_block_id_lookup(block_lengths : &[i32]) -> Vec<i32>{
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