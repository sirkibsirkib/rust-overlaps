

pub mod solutions{
    use std::hash::{Hash, SipHasher, Hasher};
    use std::collections::HashMap;

    #[derive(Hash,PartialEq, Eq, Debug)]
    pub enum Orientation{
        Normal,
        Reversed,
    }

    //NOT oriented
    #[derive(Hash,PartialEq, Eq, Debug)]
    pub struct Candidate{
        pub id_b : i32,
        pub overlap_a : i32,
        pub overlap_b : i32,
        pub overhang_right_b : i32,
    }

    //oriented
    pub struct Solution{
        pub id_a : i32,
        pub id_b : i32,
        pub orientation : Orientation,
        pub overlap_a : i32,
        pub overlap_b : i32,
        pub overhang_left_a : i32,
        pub overhang_right_b : i32,
        pub errors : u32,
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
        pub id2name : HashMap<i32, String>,
        pub id2str_in_s : HashMap<i32, Vec<u8>>,
        pub bdmap_index_id : BidirMap<i32, i32>,
    }
    impl Maps{
        pub fn num_strings(&self) -> i32{
            self.id2str_in_s.len() as i32
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

pub mod string_walk{
    pub trait Walkable<'a>{
        fn read(&self) -> u8;
        fn can_read(&self) -> bool;
        fn advance(&mut self);
        fn new(&'a [u8]) -> Self;
    }

    #[derive (Clone)]
    pub struct ForwardWalker<'a>{
        src : &'a [u8],
        next_position : usize,
    }

    #[derive (Clone)]
    pub struct BackwardWalker<'a>{
        src : &'a [u8],
        next_position : usize,
    }

    impl<'a> Walkable<'a> for ForwardWalker<'a>{
        fn new(src : &'a [u8]) -> ForwardWalker<'a>{
            ForwardWalker{src:src, next_position:0}
        }

        fn read(&self) -> u8{
            self.src[self.next_position]
        }
        fn can_read(&self) -> bool{
            self.next_position < self.src.len()
        }
        fn advance(&mut self){
            self.next_position += 1;
        }
    }

    impl<'a> Walkable<'a> for BackwardWalker<'a>{
        fn new(src : &'a [u8]) -> BackwardWalker<'a>{
            BackwardWalker{src:src, next_position:src.len()-1}
        }
        fn read(&self) -> u8{
            self.src[self.next_position]
        }
        fn can_read(&self) -> bool{
            self.next_position < self.src.len() && self.next_position > 0
        }
        fn advance(&mut self){
            self.next_position -= 1;
        }
    }
}