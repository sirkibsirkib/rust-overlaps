extern crate bio;

#[macro_use]
extern crate clap;

use std::io;
use bio::io::fasta;
use std::collections::HashMap;
use std::fs::File;

//use clap::{App, Arg, SubCommand};

#[derive(Debug)]
struct Config{
    input : String,
    output : String,
    reversals : bool,
    inclusions : bool,
    edit_distance : bool,
}

fn read_file(filename : &str) -> Result<HashMap<i32, Vec<u8>>, io::Error> {
    let mut id2name     : HashMap<i32, String> = HashMap::new();
    let mut id2str_in_s : HashMap<i32, Vec<u8>> = HashMap::new();
    let mut next_id = 0;

    let f = File::open(filename).expect("Unable to open file");
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
    Ok(id2str_in_s)
}

fn main(){
    let matches = clap_app!(ASPOPsolver =>
        (version: "1.0")
        (author: "Christopher Esterhuyse <christopher.esterhuyse@gmail.com>")
        (about: "Finds approximate suffix prefix overlaps")
        (@arg IN_PATH: +required +takes_value "Path to the input fasta file")
        (@arg OUT_PATH: +required +takes_value "Path of desired output file")
        (@arg reversals: -r "Enables reversals of S strings")
        (@arg inclusions: -i "Enables finding of inclusion solutions")
        (@arg edit_distance: -e "Uses edit distance instead of Hamming")
    ).get_matches();

    let config = Config{
        input  : matches.value_of("IN_PATH").unwrap().to_owned(),
        output : matches.value_of("OUT_PATH").unwrap().to_owned(),
        reversals : if matches.occurrences_of("reversals") >= 1 {true} else {false},
        inclusions : if matches.occurrences_of("inclusions") >= 1 {true} else {false},
        edit_distance : if matches.occurrences_of("edit_distance") >= 1 {true} else {false},
    };

    println!("{:?}", &config);

    let filename = "data/basic.fasta";
    let r = read_file(&filename);
    match r {
        Ok(_) => println!("YAY"),
        Err(_) => println!("BOO"),
    }
}