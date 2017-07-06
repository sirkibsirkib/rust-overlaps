use std::fmt::{Display, Debug};

pub mod kucherov;
pub mod valimaki;
pub type Mode = Box<IsMode>;

/*
"interface" for new filtering and partition schemes.
1. Create any struct that implements these functions
2. Add your new struct to the code in setup.rs so that the solver will use it when the arg is used
*/
pub trait IsMode: Sync + Display + Debug {

    /*
    filtering scheme. Return the number of permitted errors for a query search node with given properties
    "completed_blocks" : number of fully-matched blocks so far in THIS query search
    "patt_blocks" : number of blocks the pattern this search is for is divided into
    "blind_blocks" : number of blocks to the LEFT of this search i.e. not involved in the search
    */
    fn filter_func(&self, completed_blocks : i32, patt_blocks : i32, blind_blocks : i32) -> i32;

    // partition scheme. For a pattern of given length and alg parameters, return a vector of block lengths. order will be respected
    fn get_block_lengths(&self, patt_len : i32, err_rate : f32, thresh : i32) -> Vec<i32>;

    // return true IFF a node with the properties represented by the args should generate candidates
    fn candidate_condition(&self,generous_overlap_len : i32, completed_blocks : i32, thresh : i32, errors : i32 ) -> bool;

    // The pattern will only create query searches for pattern-block-sequence suffixes of this length or more
    fn get_fewest_suff_blocks(&self) -> i32;

    // Used by testing.rs for the cargo testing
    fn get_guaranteed_extra_blocks(&self) -> i32;
}
/*
Add your custom modes in this switch statement so that
they will be used when the solver is run with the appropriate -m flag arg.
*/
pub fn get_mode(arg : &str) -> Mode {
    let tokens : Vec<&str> = arg.split('_').collect();
    if tokens.len() == 0 {
        panic!("")
    }
    let mode_args = &tokens[1..];
    match tokens[0] {
        "valimaki" => Box::new(valimaki::ValimakiMode::new()),
        "kucherov" => Box::new(kucherov::KucherovMode::new(mode_args)),
        /*
        NEW MODE OPTIONS GO IN THIS BLOCK
        CATCH the name you want it to be associated with, whatever you like.
        return a box contining your IsMode-implementing struct like this:
            your_mod_rs_file::YourStruct::new(mode_args)
        ("IsMode" trait is defined above)
        You can also leave out the mode_args if your new() is defined as requiring no parameter.
        */






        // YOUR MODES GO HERE ^^^^
        _ => panic!("No mode with the given name found!"),
    }
}


pub fn default_mode() -> Mode {
    Box::new(kucherov::KucherovMode::new(&vec!["2"]))
}