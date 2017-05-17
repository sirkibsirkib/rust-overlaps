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

mod my_io;
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
    let config = my_io::parse_run_args();
    if config.verbose {println!("OK interpreted config args.\n Config : {:#?}", &config)};
    let (maps, text) = prepare::read_and_prepare(&config.input, &config).expect("Couldn't interpret data.");
    if config.verbose {println!("OK read and mapped fasta input.")};

    solve(&text, &config, &maps);
    if config.verbose {println!("OK done.")};
}

fn solve(text : &Vec<u8>, config : &Config, maps : &Maps){


    let alphabet = Alphabet::new(b"$ACGTN");
    let sa = suffix_array(&text);
    let bwt = bwt(&text, &sa);
    let less = less(&bwt, &alphabet);
    let occ = Occ::new(&bwt, 3, &alphabet);
    let fm = FMIndex::new(&bwt, &less, &occ);

    let pattern = b"CATTTGTA";
    println!("pattern : {}", String::from_utf8_lossy(pattern));

    let sai = fm.backward_search(pattern.iter());
    let positions = sai.occ(&sa);
    let mut positions = positions.to_owned();
    positions.sort();
    println!("interval {:?}", &sai);
    println!("positions {:?}", &positions);
    println!("text:\n{}", String::from_utf8_lossy(&text));
    for p in positions{
        let mut res = String::new();
        for _ in 0..p{
            res.push(' ');
        }
        res.push_str(&String::from_utf8_lossy(pattern));
        println!("{}", res);
    }
    let f = File::create(&config.output).expect("Unable to create output file");
    let mut wrt_buf = BufWriter::new(f);
    wrt_buf.write_all("idA\tidB\tO\tOHA\tOHB\tOLA\tOLB\tK\tCIGAR\n".as_bytes()).expect("ech");

    pipeline(
        "pipelinexyz",   // name of the pipeline for logging
         1,              // number of worker threads
         maps.id2str_in_s.keys(),   // iterator with work items

         |n| fm.generate_candidates(pattern, config, maps, *n, &sa), // computation to apply in parallel to work items
         |r| {           // aggregation to apply to work results
             println!("writing all {:?} candidates.", r.len());
             for c in r.iter() {
                 println!("WRITING CAND");
                 wrt_buf.write((format!("{}\t{}\t{}\t{}\n", c.id_b, c.overlap_a, c.overlap_b, c.overhang_right_b)).as_bytes());
             }
         }
    );
}

impl<DBWT: DerefBWT + Clone, DLess: DerefLess + Clone, DOcc: DerefOcc + Clone> GeneratesCandidates
                    for FMIndex<DBWT, DLess, DOcc> {

}