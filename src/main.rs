extern crate bidir_map;

extern crate bio;
use bio::data_structures::bwt::{DerefBWT, DerefOcc, DerefLess};
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::{FMIndex, FMIndexable};
use bio::data_structures::suffix_array::suffix_array;
use bio::data_structures::suffix_array::SuffixArray;
use bio::alphabets::Alphabet;

#[macro_use]
extern crate clap;

extern crate cue;
use cue::pipeline;

use std::fs::File;
use std::io::{Write, BufWriter};


/////////////////////////////////////

//mod algorithm_modes;
//use algorithm_modes::kucherov as filter_mode;

mod setup;
mod prepare;
mod search;
mod verification;

mod structs;
use structs::solutions::*;
use structs::run_config::*;

mod algorithm_modes;


use std::collections::HashSet;

use search::GeneratesCandidates;


pub static ALPH : &'static [u8] = b"ACGNT";
pub static READ_ERR : u8 = b'N';

fn main(){
    let config = setup::parse_run_args();
    if config.verbose {println!("OK interpreted config args.\n Config : {:#?}", &config)};
    let maps = prepare::read_and_prepare(&config.input, &config).expect("Couldn't interpret data.");
    if config.verbose {println!("OK read and mapped fasta input.")};

    solve(&config, &maps);
    println!("OK done.");
}

fn solve(config : &Config, maps : &Maps){
    let alphabet = Alphabet::new(b"$ACGTN");
    let sa = suffix_array(&maps.text);
    let bwt = bwt(&maps.text, &sa);
    let less = less(&bwt, &alphabet);
    let occ = Occ::new(&bwt, 3, &alphabet);
    let fm = FMIndex::new(&bwt, &less, &occ);


//    let test_patt = b"ABCD";

//    println!("pattern : {}", String::from_utf8_lossy(test_patt));

//    let sai = fm.backward_search(test_patt.iter());
//    let positions = sai.occ(&sa);
//    let mut positions = positions.to_owned();
//    positions.sort();
//    println!("interval {:?}", &sai);
//    println!("positions {:?}", &positions);
//    println!("text:\n{}", String::from_utf8_lossy(&maps.text));
//    for p in positions{
//        let mut res = String::new();
//        for _ in 0..p{
//            res.push(' ');
//        }
//        res.push_str(&String::from_utf8_lossy(test_patt));
//        println!("{}", res);
//    }
    let f = File::create(&config.output).expect("Unable to create output file");
    let mut wrt_buf = BufWriter::new(f);
    wrt_buf.write_all("idA\tidB\tO\t-LA\tRB-\tOLA\tOLB\tK\tCIGAR\n".as_bytes()).expect("ech");

    //don't need to search with reversed strings as patterns. would find redundant solutions
    let id_iterator = IDIterator{
        num_ids : maps.num_ids,
        next : 0,
        step : if config.reversals {2} else {1},
    };
//    for id in id_iterator {
//        let pattern = maps.get_string(id);
//        let cands = fm.generate_candidates(pattern, config, maps, id, &sa);
//        let solutions = verification::verify_all(id, cands, config, maps);
//        for sol in solutions{
//            write_solution(&mut wrt_buf, sol, maps);
//        }
//    }
    if config.sorted{
        panic!("Sorted output not implemented!");
    }
    if config.verbose{println!("OK spawning {} worker threads.", config.worker_threads);}
    pipeline(
        "pipelinexyz",              // name of the pipeline for logging
         config.worker_threads,     // number of worker threads
         id_iterator,                // iterator with work items
         |id_a| solve_an_id(config, maps, id_a, &sa, &fm), // computation to apply in parallel to work items
         |solutions| {                 // aggregation to apply to work results
             for sol in solutions {
                 write_solution(&mut wrt_buf, sol, maps);
             }
         }
    );
}

#[inline]
fn solve_an_id<DBWT: DerefBWT + Clone, DLess: DerefLess + Clone, DOcc: DerefOcc + Clone>
        (config : &Config, maps : &Maps, id_a : usize, sa : &SuffixArray, fm : &FMIndex<DBWT, DLess, DOcc>)
                -> HashSet<Solution>{
    let candidates = fm.generate_candidates(maps.get_string(id_a), config, maps, id_a, sa);
    verification::verify_all(id_a, candidates, config, maps)
}



#[inline]
fn write_solution(buf : &mut BufWriter<File>, s : Solution, maps : &Maps){
    let formatted = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                            maps.id2name_vec.get(s.id_a).expect("IDA IN THERE"),
                            maps.id2name_vec.get(s.id_b).expect("IDB IN THERE"),
                            if s.orientation==Orientation::Normal{"N"}else{"I"},
                            s.overhang_left_a,
                            s.overhang_right_b,
                            s.overlap_a,
                            s.overlap_b,
                            s.errors,
                            s.cigar,
    );
    buf.write(formatted.as_bytes()).is_ok();
}

struct IDIterator{
    num_ids : usize,
    next : usize,
    step : usize,
}

impl Iterator for IDIterator {
    type Item = usize;

    //inline means compiler will not do a proper function call :) but this becomes more like a macro
    #[inline]
    fn next(&mut self) -> Option<usize> {
        if self.next < self.num_ids{
            self.next += self.step;
            Some(self.next - self.step)
        }else{
            None
        }
    }
}

impl<DBWT: DerefBWT + Clone, DLess: DerefLess + Clone, DOcc: DerefOcc + Clone> GeneratesCandidates
                    for FMIndex<DBWT, DLess, DOcc> {

}