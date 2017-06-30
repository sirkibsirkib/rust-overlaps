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

    // Used by testing.rs for the cargo testing
    fn get_guaranteed_extra_blocks(&self) -> i32;

    // The pattern will only create query searches for pattern-block-sequence suffixes of this length or more
    fn get_fewest_suff_blocks(&self) -> i32;


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
}

