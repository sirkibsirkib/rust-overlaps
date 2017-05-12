extern crate bio;

use structs::solutions::*;
use structs::run_config::*;


use bio::alphabets;
use bio::data_structures::suffix_array::suffix_array;
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::{FMIndex, FMIndexable};
use bio::io::fastq;
use bio::alignment::distance::*;

fn step1(text : &str, config : &Config, maps : &Maps){
    //TODO thread_pool_1
    //TODO iterate over correct IDs that respect reversals
    for p_id in maps.id2str_in_s.keys(){
        let s = maps.id2str_in_s.get(p_id);
        let (block_lengths, filters) =
            mode::block_lengths_and_filters(s.len(), config.err_rate, config.thresh);
        let mut p_i_start = 0;
        for block_length in block_lengths{
            let pattern : &str = s[p_i_start..];
            //TODO put in the rust bio shet
            let (fwd_indices, bwd_indices) = search(p_id, pattern, config); //returns cand tups
            let raw_candidates = interpret_raw_candidates(fwd_indices, bwd_indices);
            write_raw_candidates(p_id, raw_candidates);
            p_i_start += block_length;
        }
    }
}

fn interpret_raw_candidates(fwd_indices : [i32], bwd_indices : [i32]) -> HashSet<RawCandidate>{
    //TODO interpret the output of rustbio. return candidate
}

fn write_raw_candidates(p_id : i32, raw_candidates : &HashSet<RawCandidate>){
    let filepath = TEMP_DIR.to_path_buf().push(p_id.to_string().push('h'));
    let home_f = File::create(filepath).expect("Unable to create home file");
    let mut wrt_buf = BufWriter::new(home_f);
    for rc in raw_candidates.into_iter() {
        write!(&mut wrt_buf, "{},{},{},{}\n", rc.id_b, rc.overlap_a, rc.overlap_b, rc.overhang_right_b);
    }
}
