extern crate bidir_map;

extern crate bio;
use bio::data_structures::bwt::{DerefBWT, DerefOcc, DerefLess};
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::FMIndex;
use bio::data_structures::suffix_array::suffix_array;
use bio::data_structures::suffix_array::RawSuffixArray;
use bio::alphabets::Alphabet;

#[macro_use]
extern crate clap;

extern crate cue;

use std::fs::File;
use std::io::{Write, BufWriter};
use std::collections::HashSet;
use std::time::Instant;

/////////////////////////////////////

mod setup;
mod prepare;
mod search;
mod verification;
mod structs;
mod algorithm_modes;

use structs::solutions::*;
use structs::run_config::*;
use search::GeneratesCandidates;


pub static ALPH : &'static [u8] = b"ACGNT";
pub static READ_ERR : u8 = b'N';

fn main(){
    let start_time = Instant::now();
    let config = setup::parse_run_args();
    if config.verbose {println!("OK interpreted config args.\n{:#?}", &config)};
    let maps = prepare::read_and_prepare(&config.input, &config).expect("Couldn't interpret data.");
    if config.verbose {println!("OK read and mapped fasta input.")};
    solve(&config, &maps);
    println!("OK done. Process finished in {:?}.", &Instant::elapsed(&start_time));
}

fn solve(config : &Config, maps : &Maps){
    let alphabet = Alphabet::new(b"$ACGTN");
    let sa = suffix_array(&maps.text);
    let bwt = bwt(&maps.text, &sa);
    let less = less(&bwt, &alphabet);
    let occ = Occ::new(&bwt, 3, &alphabet);
    let fm = FMIndex::new(&bwt, &less, &occ);
    let f = File::create(&config.output).expect("Unable to create output file");
    let mut wrt_buf = BufWriter::new(f);
    wrt_buf.write_all("idA\tidB\tO\t-LA\tRB-\tOLA\tOLB\tK\tCIGAR\n".as_bytes()).expect("ech");

    let id_iterator = IDIterator{
        num_ids : maps.num_ids,
        next : 0,
        //skip reversed input strings (would find redundant solutions)
        step : if config.reversals {2} else {1},
    };
    let mut complete_solution_list : Vec<Solution> = Vec::new(); //only used in event output sorting is desired
    if config.verbose{println!("OK spawning {} worker threads.", config.worker_threads);}

    //multithreaded part
    cue::pipeline(
        "overlap_pipeline",          // name of the pipeline for logging
         config.worker_threads,      // number of worker threads
         id_iterator,                // iterator with work items
         |id_a|                      // computation to apply in parallel to work items
            solve_an_id(config, maps, id_a, &sa, &fm),
         |solutions| {               // aggregation to apply to work results
             if config.sorted {
                 //workers ==> solutions --> sorted_solutions --> out
                 for sol in solutions {complete_solution_list.push(sol);}
             }else {
                 //workers ==> out
                 for sol in solutions {write_solution(&mut wrt_buf, &sol, maps, config);}
             }
         }
    );

    //sequential part
    if config.sorted{
        if config.verbose{println!("OK output list sorted.");}
        complete_solution_list.sort();
        for sol in complete_solution_list.iter(){
            write_solution(&mut wrt_buf, sol, maps, config);
        }
    }
    if config.verbose{println!("OK output file {} written.", config.output);}
}

#[inline]
fn solve_an_id<DBWT: DerefBWT + Clone, DLess: DerefLess + Clone, DOcc: DerefOcc + Clone>
        (config : &Config, maps : &Maps, id_a : usize, sa : &RawSuffixArray, fm : &FMIndex<DBWT, DLess, DOcc>)
                -> HashSet<Solution>{
    let candidates = fm.generate_candidates(maps.get_string(id_a), config, maps, id_a, sa);
    verification::verify_all(id_a, candidates, config, maps)
}


#[inline]
fn write_solution(buf : &mut BufWriter<File>, s : &Solution, maps : &Maps, config : &Config){
    let formatted = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                            maps.get_name_for(s.id_a),
                            maps.get_name_for(s.id_b),
                            if s.orientation==Orientation::Normal{"N"}else{"I"},
                            s.overhang_left_a,
                            s.overhang_right_b,
                            s.overlap_a,
                            s.overlap_b,
                            s.errors,
                            s.cigar,
    );
    buf.write(formatted.as_bytes()).is_ok();
    if config.print{
        let a = &String::from_utf8_lossy(maps.get_string(s.id_a));
        let b = &String::from_utf8_lossy(maps.get_string(s.id_b));
        let a_name = maps.get_name_for(s.id_a);
        let b_name = maps.get_name_for(s.id_b);
        if s.overhang_left_a > 0{
            let space = &std::iter::repeat(" ").take(s.overhang_left_a as usize).collect::<String>();
            println!(" '{}':\t{}\n [{}]:\t{}{}\n", a_name, a, b_name, space, b);
        }else{
            let space = &std::iter::repeat(" ").take((-s.overhang_left_a) as usize).collect::<String>();
            println!(" '{}':\t{}{}\n [{}]:\t{}\n", a_name, space, a, b_name, b);
        }
    }
}

#[derive(Debug, Clone)]
struct IDIterator{
    num_ids : usize,
    next : usize,
    step : usize,
}

impl Iterator for IDIterator {
    type Item = usize;

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
    // empty
}