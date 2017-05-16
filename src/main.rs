#[allow(unused_imports)]


extern crate bidir_map;
use bidir_map::BidirMap;

extern crate bio;
use bio::io::fasta;
use bio::data_structures::fmindex::*;
use bio::alphabets::dna;
use bio::data_structures::suffix_array::RawSuffixArray;
use bio::data_structures::bwt::{BWT, DerefBWT, DerefOcc, DerefLess};

#[macro_use]
extern crate clap;
extern crate cue;

use std::io;
use std::collections::HashMap;
use std::path::Path;
use std::fs;
use std::collections::HashSet;
use std::fs::File;
use std::io::{Write, BufWriter};


/////////////////////////////////////

mod algorithm_modes;
use algorithm_modes::kucherov as mode;

mod structs;
use structs::solutions::*;
use structs::run_config::*;

//mod p1;
//mod p2;


pub static ALPH : &'static [u8] = b"ACGNT";
pub static INDEX_ALPH : &'static [u8] = b"ACGNT$";
pub static READ_ERR : u8 = b'N';


fn parse_run_args() -> Config{
    let matches = clap_app!(ASPOPsolver =>
        (version: "1.0")
        (author: "Christopher Esterhuyse <christopher.esterhuyse@gmail.com>")
        (about: "Finds approximate suffix prefix overlaps")
        (@arg IN_PATH: +required +takes_value "Path to the input fasta file")
//        (@arg TMP_PATH: +required +takes_value "Path for the program to make temp files")
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
    let config = parse_run_args();
    if config.verbose {println!("OK interpreted config args.\n Config : {:#?}", &config)};
    let (maps, text) = read_and_prepare(&config.input, &config).expect("Couldn't interpret data.");
    if config.verbose {println!("OK read and mapped fasta input.")};

    solve(&text, &config, &maps);
    if config.verbose {println!("OK done.")};
}


fn read_and_prepare(filename : &str, config : &Config) -> Result<(Maps, Vec<u8>), io::Error> {
    let mut next_id = 0;
    let mut id2name : HashMap<i32, String> = HashMap::new();
    let mut id2str_in_s : HashMap<i32, Vec<u8>> = HashMap::new();
    let f = File::open(filename)
        .expect(&format!("Failed to open input file at {:?}\n", filename));
    let reader = fasta::Reader::new(f);
    for record in reader.records() {
        let record = record?;
        if let Some(name) = record.id(){
            print!("ID:{:?}, seq {:?}", name, String::from_utf8_lossy(record.seq()));
            println!();
            id2name.insert(next_id, name.to_owned());
            id2str_in_s.insert(next_id, record.seq().to_vec());
            next_id += 1;
        }
    }

    let (text, bdmap_index_id) = make_text(config.reversals, &id2str_in_s);

    id2name.shrink_to_fit();
    id2str_in_s.shrink_to_fit();

    let maps = Maps{
        id2name : id2name,
        id2str_in_s : id2str_in_s,
        bdmap_index_id : bdmap_index_id,
    };

    Ok((maps, text))
}

fn make_text(reverse : bool, strings : &HashMap<i32, Vec<u8>>)
            -> (Vec<u8>, BidirMap<i32, i32>){
    //TODO need to consider reversals
    let mut bdmap_index_id : BidirMap<i32, i32> = BidirMap::new();
    let mut end_dollar2id : HashMap<i32, i32> = HashMap::new();
    let mut text : Vec<u8> = Vec::new();
    for id in 0..(strings.len() as i32){
        let index = text.len() as i32;
        bdmap_index_id.insert(index, id);
        let s = strings.get(&id).expect("Can't find string with that ID!");
        //TODO reverse string + flip symbols sometimes
        text.push('$' as u8);
        text.extend(s);
    }
    text.push('$' as u8);
    text.shrink_to_fit();
    (text, bdmap_index_id)
}

fn solve(text : &Vec<u8>, config : &Config, maps : &Maps){
    use cue::pipeline;
    use bio::data_structures::bwt::{bwt, less, Occ};
    use bio::data_structures::fmindex::{FMIndex, FMIndexable};
    use bio::data_structures::suffix_array::suffix_array;
    use bio::alphabets::dna;

    let alphabet = dna::n_alphabet();
    let sa = suffix_array(&text);
    let bwt = bwt(&text, &sa);
    let less = less(&bwt, &alphabet);
    let occ = Occ::new(&bwt, 3, &alphabet);
    let fm = FMIndex::new(&bwt, &less, &occ);

    let pattern = b"CAT";
    println!("pattern : {}", String::from_utf8_lossy(pattern));

    let sai = fm.backward_search(pattern.iter());
    let positions = sai.occ(&sa);
    let mut positions = positions.to_owned();
    positions.sort();
    println!("{:?}", &positions);
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
         |n| work(*n, maps.id2str_in_s.get(n).expect("WTF"), config, &fm, maps, &sa), // computation to apply in parallel to work items
         |r| {           // aggregation to apply to work results
             for c in r.iter() {
                 wrt_buf.write((format!("{}\t{}\t{}\t{}\n", c.id_b, c.overlap_a, c.overlap_b, c.overhang_right_b)).as_bytes());
             }
         }
    );
}

use bio::data_structures::fmindex::{FMIndexable, FMIndex};
fn work<DBWT: DerefBWT + Clone,
    DLess: DerefLess + Clone,
    DOcc: DerefOcc + Clone>(
        p_id : i32, pattern : &[u8],
        config : &Config,
        fm : &FMIndex<DBWT, DLess, DOcc>,
        maps : &Maps,
        sa : &RawSuffixArray,
        )
        -> HashSet<Candidate>{

    println!("New work with p_id {}, pattern {}", p_id, String::from_utf8_lossy(pattern));
    let patt_len = pattern.len() as i32;
    let block_lengths = alg_mode::get_block_lengths(patt_len, config.err_rate, config.thresh);
    println!("block lengtsh {:?} patt len {}", &block_lengths, patt_len);
    fm.generate_candidates(
        pattern,
        config,
        maps,
        p_id,
        &block_lengths,
        sa,
    )
}



impl<DBWT: DerefBWT + Clone, DLess: DerefLess + Clone, DOcc: DerefOcc + Clone> GeneratesCandidates
                    for FMIndex<DBWT, DLess, DOcc> {

}

use algorithm_modes::kucherov as alg_mode;
trait GeneratesCandidates : FMIndexable {

    fn generate_candidates(&self,
            pattern : &[u8],
            config : &Config,
            maps : &Maps,
            id_a : i32,
            block_lengths : &Vec<i32>,
            sa : &RawSuffixArray,
            )
            -> HashSet<Candidate> {

        let patt_len = pattern.len() as i32;
        let mut candidate_set: HashSet<Candidate> = HashSet::new();
        let max_b_len =
            if config.reversals {patt_len} else {(patt_len as f32 / (1.0 - config.err_rate)).floor() as i32};
        let max_errs = (max_b_len as f32 * config.err_rate).floor() as i32;
        let block_id_lookup = get_block_id_lookup(max_b_len, config, block_lengths);
        println!("block_id_lookup {:?}", &block_id_lookup);


        let full_interval = Interval {
            lower: 0,
            upper: self.bwt().len() -1,
        };

        let mut p_i : i32 = 0;

        // FOR EACH FILTER, from smallest prefix to total prefix
        // PREFIX FILTERS
        println!("\n\nGenerate candidates!");
        let patt_blocks : i32 = block_lengths.len() as i32;
        for (block_id, block_len) in block_lengths.iter().enumerate() {
            println!("\n\n\nnew suff of len {} for {}", block_id+1, String::from_utf8_lossy(&pattern));
            let search_constants = SearchConstants{
                pattern: pattern,
                config : config,
                maps : maps,
                block_id_lookup : &block_id_lookup,
                sa : sa,
                id_a : id_a,
                first_block_id : block_id as i32,
                patt_blocks : patt_blocks,
            };
            println!(">> id_a {} first_block_id {}", search_constants.id_a, search_constants.first_block_id);
            self.recurse_candidates(&mut candidate_set, &search_constants, 0, p_i, 0, 0, 0, &full_interval, &String::new());
            p_i += (*block_len) as i32;
        }
        candidate_set
    }

    fn recurse_candidates(&self,
                          cand_set : &mut HashSet<Candidate>,
                          cns : &SearchConstants,
                          errors : i32,
                          p_i : i32,
                          indel_balance : i32,
                          a_match_len : i32,
                          b_match_len : i32,
                          matches : &Interval,
                          debug : &str){
        println!("'{}'   recurse p_i=={}    {:?}", debug, p_i, &matches);

        if matches.lower > matches.upper{
            // range is (lower, upper)  ie. both inclusive
            return
        }

        //block_id_lookup has 1 char of padding on the end so we can avoid checking bounds
        let completed_blocks : i32 = cns.block_id_lookup[p_i as usize]- cns.first_block_id;
        let permitted_error : i32 =
            alg_mode::filter_func(completed_blocks, cns.patt_blocks);

//        let has_spare_error : bool = errors < permitted_error;

        // ADD SUFFIX CANDIDATES
        let generous_match_len = std::cmp::max(a_match_len, b_match_len);

        let cand_condition_satisfied =
            alg_mode::candidate_condition(generous_match_len, completed_blocks, cns.config.thresh, errors);

        if cand_condition_satisfied {
            println!("  !!!!!! CAND COND SATISFIED! :D {}", debug);
            let a = '$' as u8;
            let less = self.less(a);
            let dollar_interval = Interval {
                lower : less + if matches.lower > 0 { self.occ(matches.lower - 1, a) } else { 0 },
                upper : less + self.occ(matches.upper, a) - 1,
            };
            println!("interval {:?}    ==$==> {:?}", &matches, &dollar_interval);
            let positions = dollar_interval.occ(cns.sa);
            println!("positions {:?}", &positions);
            for p in positions {
                let id_b = *cns.maps.bdmap_index_id.get_by_first(&(p as i32 + 1)).expect("DOLLAR MAP BAD");
                if id_b == cns.id_a{
                    //skip obvious self-matches
                    continue
                }
//                let c = Candidate {
//                    id_b: id_b,
//                    overlap_a: a_match_len,
//                    overlap_b: b_match_len,
//                    overhang_right_b: 0,
//                };
                println!("  !!!!!! ~~~  adding FOUND CANDIDATE AT {} with {}", p, debug);
//                cand_set.insert(c);
            }
        }

        if p_i >= cns.pattern.len() as i32{
            //END OF PATTERN
            //INCLUSIONS GO HERE
            return
        }


        // walking + substitution match
        for &a in ALPH.iter(){
            let less = self.less(a);
            let next_interval = Interval{
                lower : less + if matches.lower > 0 { self.occ(matches.lower - 1, a) } else { 0 },
                upper : less + self.occ(matches.upper, a) - 1,
            };
            let p_char_real_index = (cns.pattern.len() - (p_i) as usize -1);
            let p_char = cns.pattern[p_char_real_index];
            let recurse_error =  if p_char == a && a != READ_ERR {errors} else {errors + 1};
            if recurse_error <= permitted_error{
                let mut next_debug = format!("{}{}", a as char, debug);
                self.recurse_candidates(cand_set,
                                        cns,
                                        recurse_error,
                                        p_i+1,
                                        0,
                                        a_match_len+1,
                                        b_match_len+1,
                                        &next_interval,
                                        &next_debug);
            }
        }
    }
}

#[derive(Debug)]
struct SearchConstants<'a>{
    config : &'a Config,
    maps : &'a Maps,

    block_id_lookup : &'a Vec<i32>,
    sa : &'a RawSuffixArray,
    pattern: &'a [u8],
    id_a : i32,

    first_block_id : i32,
    patt_blocks : i32,
}

fn get_block_id_lookup(max_b_len : i32, config : &Config, block_lengths : &[i32]) -> Vec<i32>{
    let mut lookup : Vec<i32> = Vec::new();
    for (id, block_length) in block_lengths.iter().enumerate() {
        for _ in 0..*block_length{
            lookup.push(id as i32);
        }
    }
    let last_index = block_lengths.len() as i32 - 1;
    while (lookup.len() as i32) < max_b_len + 1{
        lookup.push(last_index);
    }
    lookup.shrink_to_fit();
    lookup
}


//    fn garbage<'b, P: Iterator<Item = &'b u8> + DoubleEndedIterator>(&self,
//                                                                             pattern: P)
//                                                                             -> Interval {
//        let (mut l, mut r) = (0, self.bwt().len() - 1);
//        for &a in pattern.rev() {
//            let less = self.less(a);
//            l = less + if l > 0 { self.occ(l - 1, a) } else { 0 };
//            r = less + self.occ(r, a) - 1;
//        }
//
//        Interval {
//            lower: l,
//            upper: r + 1,
//        }
//    }
//}