#[allow(unused_imports)]


#[macro_use]
extern crate clap;

use std::io;
use std::collections::HashMap;
use std::path::Path;
use std::fs;
use std::collections::HashSet;
use std::fs::File;
use std::io::{Write, BufWriter};

extern crate bio;
use bio::io::fasta;


extern crate bidir_map;
use bidir_map::BidirMap;

/////////////////////////////////////

mod algorithm_modes;
use algorithm_modes::kucherov as mode;

mod structs;
use structs::solutions::*;
use structs::run_config::*;

//mod p1;
//mod p2;

#[allow(dead_code)]


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
    let (maps, text) = read_and_prepare(&config.input, &config).expect("Couldn't interpret data");
    if config.verbose {println!("OK read and mapped fasta input")};
    if !Path::new(structs::TEMP_DIR).exists(){
        let r = fs::create_dir(structs::TEMP_DIR);
        assert!(r.is_ok());
        println!("Created /tmp/.");
    } else{
        println!("/tmp/ already here");

    }

    println!("{:?}", String::from_utf8_lossy(&text));

//    p1::step1(text, config, maps);
//    let all_solutions =
//        p2::step2(text, config, maps);
//    write_solutions(&all_solutions);
}


fn read_and_prepare(filename : &str, config : &Config) -> Result<(Maps, Vec<u8>), io::Error> {
    let mut next_id = 0;
    let mut id2name : HashMap<usize, String> = HashMap::new();
    let mut id2str_in_s : HashMap<usize, Vec<u8>> = HashMap::new();
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

fn make_text(reverse : bool, strings : &HashMap<usize, Vec<u8>>)
            -> (Vec<u8>, BidirMap<usize, usize>){
    //TODO need to consider reversals
    let mut bdmap_index_id : BidirMap<usize, usize> = BidirMap::new();
    let mut text : Vec<u8> = Vec::new();
    for id in 0..strings.len(){
        let index = text.len();
        bdmap_index_id.insert(index, id);
        let s = strings.get(&id).expect("Can't find string with that ID!");
        //TODO reverse string + flip symbols sometimes
        text.extend(s);
        text.push('$' as u8);
    }
    text.shrink_to_fit();
    (text, bdmap_index_id)
}

fn write_solutions(solution_sets : &Vec<HashSet<Solution>>, config : &Config){
    let f = File::create(&config.output).expect("Unable to create output file");
    let mut wrt_buf = BufWriter::new(f);
    wrt_buf.write_all("idA\tidB\tO\tOHA\tOHB\tOLA\tOLB\tK\tCIGAR\n".as_bytes())
        .expect("Couldn't write first output line");
    for solution_set in solution_sets.iter() {
        for sol in solution_set.iter() {
            let s = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                            sol.id_a,
                            sol.id_b,
                            sol.orientation,
                            sol.overhang_left_a,
                            sol.overhang_right_b,
                            sol.overlap_a,
                            sol.overlap_b,
                            sol.errors,
                            sol.cigar,
            );
            //TODO write! not working?
            wrt_buf.write(s.as_bytes())
                .expect("Failed to write solution");
        }
    }
}