extern crate bidir_map;

extern crate bio;
use bio::data_structures::bwt::{DerefBWT, DerefOcc, DerefLess};
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::{FMIndex, FMIndexable};
use bio::data_structures::suffix_array::suffix_array;
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
use algorithm_modes::kucherov::get_block_lengths;


use search::GeneratesCandidates;


pub static ALPH : &'static [u8] = b"ACGNT";
pub static READ_ERR : u8 = b'N';

fn main(){
    let config = setup::parse_run_args();
    if config.verbose {println!("OK interpreted config args.\n Config : {:#?}", &config)};
    let maps = prepare::read_and_prepare(&config.input, &config).expect("Couldn't interpret data.");
    if config.verbose {println!("OK read and mapped fasta input.")};

    solve(&config, &maps);
    if config.verbose {println!("OK done.")};
}

fn solve(config : &Config, maps : &Maps){
    let alphabet = Alphabet::new(b"$ACGTN");
    let sa = suffix_array(&maps.text);
    let bwt = bwt(&maps.text, &sa);
    let less = less(&bwt, &alphabet);
    let occ = Occ::new(&bwt, 3, &alphabet);
    let fm = FMIndex::new(&bwt, &less, &occ);


    let test_patt = b"ABCD";

    println!("pattern : {}", String::from_utf8_lossy(test_patt));

    let sai = fm.backward_search(test_patt.iter());
    let positions = sai.occ(&sa);
    let mut positions = positions.to_owned();
    positions.sort();
    println!("interval {:?}", &sai);
    println!("positions {:?}", &positions);
    println!("text:\n{}", String::from_utf8_lossy(&maps.text));
    for p in positions{
        let mut res = String::new();
        for _ in 0..p{
            res.push(' ');
        }
        res.push_str(&String::from_utf8_lossy(test_patt));
        println!("{}", res);
    }
    let f = File::create(&config.output).expect("Unable to create output file");
    let mut wrt_buf = BufWriter::new(f);
    wrt_buf.write_all("idA\tidB\tO\tOHA\tOHB\tOLA\tOLB\tK\tCIGAR\n".as_bytes()).expect("ech");

    pipeline(
        "pipelinexyz",   // name of the pipeline for logging
         1,              // number of worker threads
         0..maps.num_ids,   // iterator with work items

         |n| (fm.generate_candidates(maps.get_string(n), config, maps, n, &sa), n), // computation to apply in parallel to work items
         |(r, n)| {           // aggregation to apply to work results
             println!("writing all {:?} candidates.", r.len());
             let solutions = verification::verify_all(n, r, config, maps);
             for sol in solutions {
                 println!("WRITING CAND");
                 let formatted = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                                 sol.id_a,
                                 sol.id_b,
                                 if sol.orientation==Orientation::Normal{"N"}else{"I"},
                                 sol.overlap_a,
                                 sol.overlap_b,
                                 sol.overhang_left_a,
                                 sol.overhang_right_b,
                                 sol.errors,
                                 sol.cigar,
                 );
                 wrt_buf.write(formatted.as_bytes());
             }
         }
    );
}

impl<DBWT: DerefBWT + Clone, DLess: DerefLess + Clone, DOcc: DerefOcc + Clone> GeneratesCandidates
                    for FMIndex<DBWT, DLess, DOcc> {

}