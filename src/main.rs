#[macro_use]
extern crate clap;

use std::io;
use bio::io::fasta;
use std::collections::HashMap;
use std::path::Path;
use std::time::Duration;
use std::thread;
use std::fs::File;
use std::io::Write;
use std::sync::Mutex;
use bio::alphabets;
use bio::data_structures::suffix_array::suffix_array;
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::{FMIndex, FMIndexable};
use bio::io::fastq;
use bio::alignment::distance::*;

use bidir_map::BidirMap;

/////////////////////////////////////

mod algorithm_modes;
use algorithm_modes::kucherov as mode;

mod structs;
use structs::solutions::*;
use structs::run_config::*;

mod p1;
mod p2;

const TEMP_DIR : fs::Path = Path::new("/tmp/");

fn read_file(maps : &mut Maps, filename : &str) -> Result<(), io::Error> {
    let mut next_id = 0;

    let f = File::open(filename)
        .expect(&format!("Failed to open input file at {:?}\n", filename));
    let reader = fasta::Reader::new(f);
    for record in reader.records() {
        let record = record?;
        if let Some(name) = record.id(){
//            print!("ID:{:?}, seq {:?}", name, String::from_utf8_lossy(record.seq()));
//            println!();
            maps.id2name
                .insert(next_id, name.to_owned());
            maps.id2str_in_s
                .insert(next_id, record.seq().to_vec());
            next_id += 1;
        }
    }
    maps.id2name.shrink_to_fit();
    maps.id2str_in_s.shrink_to_fit();
    Ok(())
}


fn parse_run_args() -> Config{
    let matches = clap_app!(ASPOPsolver =>
        (version: "1.0")
        (author: "Christopher Esterhuyse <christopher.esterhuyse@gmail.com>")
        (about: "Finds approximate suffix prefix overlaps")
        (@arg IN_PATH: +required +takes_value "Path to the input fasta file")
        (@arg OUT_PATH: +required +takes_value "Path of desired output file")
        (@arg ERR_RATE: +required +takes_value "The max rate of errors in an overlap")
        (@arg THRESH: +required +takes_value "Shortest allowed length of an overlap")
        (@arg WORKER_THREADS: +required +takes_value "Number of worker threads used")
        (@arg reversals: -r "Enables reversals of S strings")
        (@arg inclusions: -i "Enables finding of inclusion solutions")
        (@arg edit_distance: -e "Uses edit distance instead of Hamming")
        (@arg verbose: -v "Prints completed steps of the run process")
    ).get_matches();

    Config{
        input  : matches.value_of("IN_PATH").unwrap().to_owned(),
        output : matches.value_of("OUT_PATH").unwrap().to_owned(),
        err_rate : matches.value_of("ERR_RATE").unwrap().parse().unwrap(),
        thresh : matches.value_of("THRESH").unwrap().parse().unwrap(),
        worker_threads: matches.value_of("WORKER_THREADS").unwrap().parse().unwrap(),
        reversals :     if matches.occurrences_of("reversals")     >= 1 {true} else {false},
        inclusions :    if matches.occurrences_of("inclusions")    >= 1 {true} else {false},
        edit_distance : if matches.occurrences_of("edit_distance") >= 1 {true} else {false},
        verbose :       if matches.occurrences_of("verbose")       >= 1 {true} else {false},
    }
}

fn main(){
    let mut maps = Maps{
        id2name : HashMap::new(),
        id2str_in_s : HashMap::new(),
        bdmap_index_id : BidirMap::new(),
    };
    let config = parse_run_args();
    if config.verbose {println!("OK interpreted config args.\n Config : {:#?}", &config)};
    read_file(&mut maps, &config.input).expect("Failed to populate string dict from file");
    if config.verbose {println!("OK read and mapped fasta input")};

    let text : String = make_text(config.reversals, maps.id2str_in_s.values(), &mut bdmap_index_id);

    if !Path::new(TEMP_DIR).exists(){
        fs::create_dir(TEMP_DIR)?;
    }

    p1::step1(text, config, maps);
    p2::step2(text, config, maps);
}


fn make_text(reverse : bool, strings : &HashMap<i32, &str>, bdmap_index_id : &mut HashMap<i32, i32>)
            -> String{
    //TODO need to consider reversals
    let mut result = String::new();
    for id in 0..strings.len(){
        let index = result.len();
        bdmap_index_id.insert(index, id);
        let s = strings.get(id);
        //TODO reverse string + flip symbols sometimes
        result.push_str(s);
        result.push('$');
    }
    result
}
