extern crate bio;

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


use bio::alphabets;
use bio::data_structures::suffix_array::suffix_array;
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::{FMIndex, FMIndexable};
use bio::io::fastq;
use bio::alignment::distance::*;

mod algorithm_modes;
use algorithm_modes::kucherov as mode;

mod sol;
use sol::solutions;
use sol::solutions::{RawCandidate, Candidate, Solution};

use bidir_map::BidirMap;

const TEMP_DIR : fs::Path = Path::new("/tmp/");


#[derive(Debug)]
struct Maps{
    id2name : HashMap<i32, String>,
    id2str_in_s : HashMap<i32, Vec<u8>>,
    bdmap_index_id : BidirMap<i32, i32>,
}
impl Maps{
    fn num_strings(&self) -> usize{
        self.id2str_in_s.len()
    }
}

#[derive(Debug)]
struct Config{
    //required
    input : String,
    output : String,
    err_rate : f32,
    thresh : i32,
    worker_threads: i32,

    //optional
    reversals : bool,
    inclusions : bool,
    edit_distance : bool,
    verbose : bool,
}


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
    step1(text, config, maps);
    step2(text, config, maps);
}

fn make_text(reverse : bool, strings : &HashMap<i32, &str>, bdmap_index_id : &mut HashMap<i32, i32>)
            -> String{
    //TODO need to consider reversals
    let mut result = String::new();
    for id in 0..strings.len(){
        let index = result.len();
        bdmap_index_id.insert(index, id);
        let s = strings.get(id);
        result.push_str(s);
        result.push('$');
    }
    result
}

fn step2(text : &str, config : &Config, maps : &Maps){
    //todo thread_pool_2
}

fn step1(text : &str, config : &Config, maps : &Maps){
    //todo thread_pool_1
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
    //TODO interpret the output of rustbio. return rawcandidate
}

fn write_raw_candidates(p_id : i32, raw_candidates : &HashSet<RawCandidate>){
    let filepath = TEMP_DIR.to_path_buf().push(p_id.to_string());
    let mut wtr = csv::Writer::from_file(Path::new(filepath));
    for rc in raw_candidates.into_iter() {
        let tup = (
            p_id,
            rc.id_b,
            rc.overlap_a,
            rc.overlap_b,
            rc.overhang_right_b,
        );
        let result = wtr.encode(tup);
        assert!(result.is_ok());
    }
    //flush writer? idk
}

fn verify_all(candidates : &HashSet<Candidate>) -> HashSet<Solution>{
    let solution_set : HashSet<Solution> = HashSet:new();
    for c in candidates{
        //TODO
        match verify()
    }
}

// incoming cands are already oriented
fn verify(id_a : i32, c : &Candidate, config : &Config, maps : &Maps) -> Option(Solution){
    let a_len = maps.id2str_in_s.get(id_a).len();
    let b_len = maps.id2str_in_s.get(c.id_b).len();


    let a_start = max(0, c.overhang_left_a);
    let a_end = a_start + c.overlap_a;
    let a_end2 = min(a_len, a_len-c.overhang_right_b);
    assert_eq!(a_end, a_end2);          //DEBUG! can remove a_end2 if this checks out

    let b_start = max(0, -c.overhang_left_a);
    let b_end = b_start + c.overlap_b;
    let b_end2 = min(b_len, b_len+c.overhang_right_b);
    assert_eq!(b_end, b_end2);          //DEBUG! can remove b_end2 if this checks out

    let a_part : &str = maps.id2str_in_s.get(id_a)  [a_start..a_end];
    let b_part : &str = maps.id2str_in_s.get(c.id_b)[b_start..b_end];

    let k : i32 = if config.edit_distance{
        hamming(a_part, b_part).expect("Can't apply hamming to mismatching lengths!")
    }else{
        levenshtein(a_part, b_part)
    };
    let k_limit : f32 = config.err_rate*(max(c.overlap_a, c.overlap_b) as f32);
    if k <= k_limit{
        Some(Solution{candidate:c, errors:k, cigar:""})
    }else{
        None()
    }
}