
extern crate clap;

use clap::{App, SubCommand};


pub mod args{

    pub struct Config{
        out_path : String,
        in_path : String,
    }

    pub fn new_config() -> Config{



        Config{out_path : "OUTP", in_path : "INP... BOYATSINYESAAAHHhhHHhh"}
    }
}