
pub static TEMP_DIR : &'static str = "/tmp/";

pub mod solutions{

    use std::hash::{Hash, SipHasher, Hasher};
    use std::collections::HashMap;

    #[derive(Hash)]
    pub enum Orientation{
        Normal,
        Inverted,
    }

    //NOT oriented
    #[derive(Hash)]
    pub struct Candidate{
        id_b : i32,
        overlap_a : i32,
        overlap_b : i32,
        overhang_right_b : i32,
    }

    //oriented
    pub struct Solution{
        id_b : i32,
        orientation : bool,
        overlap_a : i32,
        overlap_b : i32,
        overhang_left_a : i32,
        overhang_right_b : i32,
        errors : u8,
        cigar : String,
    }

    impl Hash for Solution {
        fn hash<H: Hasher>(&self, state: &mut H) {
            self.id_b.hash(state);
            self.orientation.hash(state);
            self.overlap_a.hash(state);
            self.overlap_b.hash(state);
            self.overhang_left_a.hash(state);
            self.overhang_right_b.hash(state);
            self.errors.hash(state);
            // cigar string not involved in the hash
        }
    }
}

pub mod run_config{
    extern crate bidir_map;

    use std::collections::HashMap;
    use bidir_map::BidirMap;

    #[derive(Debug)]
    pub struct Maps{
        pub id2name : HashMap<usize, String>,
        pub id2str_in_s : HashMap<usize, Vec<u8>>,
        pub bdmap_index_id : BidirMap<usize, usize>,
    }
    impl Maps{
        pub fn num_strings(&self) -> usize{
            self.id2str_in_s.len()
        }
    }

    #[derive(Debug)]
    pub struct Config{
        //required
        pub input : String,
        pub output : String,
        pub err_rate : f32,
        pub thresh : i32,
        pub worker_threads: i32,

        //optional
        pub reversals : bool,
        pub inclusions : bool,
        pub edit_distance : bool,
        pub verbose : bool,
    }
}