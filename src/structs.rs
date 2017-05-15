
pub static TEMP_DIR : &'static str = "tmp/";

pub mod solutions{

    use std::hash::{Hash, SipHasher, Hasher};
    use std::collections::HashMap;

    #[derive(Hash,PartialEq, Eq)]
    pub enum Orientation{
        Normal,
        Reversed,
    }

    //NOT oriented
    #[derive(Hash,PartialEq, Eq)]
    pub struct Candidate{
        pub id_b : usize,
        pub overlap_a : i32,
        pub overlap_b : i32,
        pub overhang_right_b : i32,
    }

    //oriented
    pub struct Solution{
        pub id_a : usize,
        pub id_b : usize,
        pub orientation : bool,
        pub overlap_a : i32,
        pub overlap_b : i32,
        pub overhang_left_a : i32,
        pub overhang_right_b : i32,
        pub errors : u8,
        pub cigar : String,
    }

    impl PartialEq for Solution {
        fn eq(&self, other: &Solution) -> bool {
            self.id_a == other.id_a
            && self.id_b == other.id_b
            && self.orientation == other.orientation
            && self.overlap_a == other.overlap_a
            && self.overlap_b == other.overlap_b
            && self.overhang_left_a == other.overhang_left_a
            && self.overhang_right_b == other.overhang_right_b
        }
    }

    // TODO why does it need to be hashable AND eq?
    impl Eq for Solution {}

    impl Hash for Solution {
        fn hash<H: Hasher>(&self, state: &mut H) {
            self.id_a.hash(state);
            self.id_b.hash(state);
            self.orientation.hash(state);
            self.overlap_a.hash(state);
            self.overlap_b.hash(state);
            self.overhang_left_a.hash(state);
            self.overhang_right_b.hash(state);
            // ERRORS and CIGAR not need not be the same
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