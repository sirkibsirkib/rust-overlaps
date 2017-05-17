use std::io;
use std::collections::HashMap;
use std::path::Path;
use std::fs;
use std::collections::HashSet;
use std::fs::File;
use std::io::{Write, BufWriter};


use bio::io::fasta;
use bidir_map::BidirMap;

use structs::run_config::*;


pub fn read_and_prepare(filename : &str, config : &Config) -> Result<(Maps, Vec<u8>), io::Error> {
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
    println!("DOLLAR MAP {:#?}", &maps);

    Ok((maps, text))
}

fn make_text(reverse : bool, strings : &HashMap<i32, Vec<u8>>)
             -> (Vec<u8>, BidirMap<i32, i32>){
    //TODO need to consider reversals
    let mut bdmap_index_id : BidirMap<i32, i32> = BidirMap::new();
    let mut end_dollar2id : HashMap<i32, i32> = HashMap::new();
    let mut text : Vec<u8> = Vec::new();
    for id in 0..(strings.len() as i32){
        text.push('$' as u8);
        let index = text.len() as i32;
        bdmap_index_id.insert(index, id);
        let s = strings.get(&id).expect("Can't find string with that ID!");
        //TODO reverse string + flip symbols sometimes
        text.extend(s);
    }
    text.push('$' as u8);
    text.push('#' as u8);
    text.shrink_to_fit();
    (text, bdmap_index_id)
}