extern crate bio;

use structs::solutions::*;
use structs::run_config::*;


use bio::alphabets;
use bio::data_structures::suffix_array::suffix_array;
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::{FMIndex, FMIndexable};
use bio::io::fastq;
use bio::alignment::distance::*;

use std::fs::File;
use std::collections::HashSet;



use algorithm_modes::kucherov as alg_mode;

use structs::TEMP_DIR;


fn step1(text : &Vec<u8>, config : &Config, maps : &Maps){
    //TODO iterate over correct IDs that respect reversals
    //TODO note that block_id_lookup is relative to entire PATTERN. the individual steps must use suff_i to offset

    for (p_id, pattern) in maps.id2str_in_s.enumerate(){
        //TODO thread_pool_1

        let patt_len = pattern.len();
        let block_lengths = alg_mode::get_block_lengths(patt_len, config.err_rate, config.thresh);
        let filter_func : &Fn(i32,i32) -> i32 = alg_mode::filter_func;
        let mut p_i_start = 0;
        let max_b_len =
            if config.reversals {patt_len} else {(patt_len as f32 / (1.0 - config.err_rate)).floor() as i32};
        let max_errs = (max_b_len as f32 * config.err_rate).floor() as i32;
        let block_id_lookup = get_block_id_lookup(max_b_len, config, &block_lengths);
        for (suff_i, block_length) in block_lengths.enumerate() {

            //TODO put in the FMIndex step

//            let cands =
//                search(p_id, suff_i, &pattern, p_i_start, config, block_id_lookup);

//            write_raw_candidates(config, p_id, &cands);
            p_i_start += block_length;
        }
    }
}


/*
consider blocks = [3|2|2]
step_err_arr  = [0|0|0|1|1|2|2]
*Note will also be extra long as needed by max_b_len to prevent out-of-bounds err
//TODO better to do out of bounds err or just lengthen the array?
*/
fn get_block_id_lookup(max_b_len : i32, config : &Config, block_lengths : &[i32]) -> Vec<u8>{
    let mut lookup : Vec<u8> = Vec::new();
    for (id, block_length) in block_lengths.enumerate() {
        for _ in 0..block_length{
            lookup.push(id);
        }
    }
    let last_index = block_lengths.len() - 1;
    while lookup.len() < max_b_len{
        lookup.push(last_index);
    }
    lookup.shrink_to_fit();
    lookup
}

fn write_raw_candidates(config : &Config, p_id : i32, cands : &HashSet<Candidate>){
    use std::io::BufWriter;

    let filepath = TEMP_DIR.to_path_buf().push(p_id.to_string());
    let f = File::create(filepath).expect("Unable to create candidate tmp file");
    let mut file_buffer = BufWriter::new(f);
    if config.reversals{
        let partner_p_id = if p_id%2==0 {p_id+1} else {p_id-1};
        for rc in cands.into_iter() {
            // firstly, dont write illogical candidates
            // secondly, don't write cands that WILL be found by a reverse
            if rc.id_b != p_id && pc.id_b != partner_p_id && (p_id < partner_p_id || rc.overhang_right_b != 0){
                write_candidate(&mut file_buffer, rc);
            }
        }
    }else{
        for rc in cands.into_iter() {
            if rc.id_b != p_id {
                write_candidate(&mut file_buffer, rc);
            }
        }
    }
}

fn write_candidate(file_buffer : &mut BufWriter, candidate : &Candidate){
    write!(&mut file_buffer, "{},{},{},{}\n",
           candidate.id_b,
           candidate.overlap_a,
           candidate.overlap_b,
           candidate.overhang_right_b
    );
}
