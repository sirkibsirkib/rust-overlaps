use bio::io::fasta;
use bidir_map::BidirMap;

use std::io;
use std::fs::File;

/////////////////////////////

use structs::run_config::{Config, Maps};

/*
builds the maps data structure from a fasta file + config
the "maps" contains most of the constant information for the run
> mappings from internal to external representations of strings ie: id-->"name"
> mappings between internal represenations ie: id<-->index(in text)
> important data ie: text
> some convenient functions ie: get &str (in the text)
*/
pub fn read_and_prepare(filename : &str, config : &Config) -> Result<(Maps), io::Error> {
    let mut text : Vec<u8> = Vec::new();
    let mut id2name_vec : Vec<String> = Vec::new();
    let mut id2index_bdmap : BidirMap<usize, usize> = BidirMap::new();

    let f = File::open(filename)
        .expect(&format!("Failed to open input file at {:?}\n", filename));
    let reader = fasta::Reader::new(f);
    let mut n_symbols_removed = 0;
    for record in reader.records() {
        let record = record?;
        if let Some(name) = record.id(){
            let id = id2name_vec.len();
            let name = name.to_owned();
            let mut str_vec = record.seq().to_vec();
            if !config.n_alphabet{
                let before_len = str_vec.len();
                str_vec.retain(|c|*c != ('N' as u8));
                if str_vec.len() < before_len{
                    n_symbols_removed += before_len - str_vec.len();
                }
            }
            str_vec.reverse();
            text.push('$' as u8);
            let index = text.len();
            id2index_bdmap.insert(id, index);
            text.extend(str_vec.clone());
            id2name_vec.push(name.clone());

            if config.reversals{
                let id = id2name_vec.len();
                str_vec.reverse();
                for i in 0..str_vec.len(){
                    str_vec[i] = complement_u8(str_vec[i]);
                }
                text.push('$' as u8);
                let index = text.len();
                id2index_bdmap.insert(id, index);
                text.extend(str_vec);
                id2name_vec.push(name);
            }
        }
    }


    if n_symbols_removed > 0 {
        println!("    WARNING\n\tOmitted {} N symbols found in input data.\n\t\
        Run without flag --no_n to use these N strings intact.", n_symbols_removed);
    }

    text.push('#' as u8);
    text.shrink_to_fit();
    id2name_vec.shrink_to_fit();

    let mut indexes : Vec<usize> = id2index_bdmap.second_col().map(|x| *x).collect();
    indexes.sort();
    indexes.shrink_to_fit();


    let maps = Maps{
        text : text,
        id2name_vec : id2name_vec,
        id2index_bdmap : id2index_bdmap,
        indexes : indexes,
    };
    Ok(maps)
}

fn complement_u8(x : u8) -> u8 {
    match x{
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'N' => b'N',
        _ => panic!("Bad string char '{}'", x as char),
    }
}
