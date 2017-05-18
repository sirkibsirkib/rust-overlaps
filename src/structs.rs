

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
        pub id_b : usize,
        pub overlap_a : usize,
        pub overlap_b : usize,
        pub overhang_right_b : i32,
        pub debug_str : String,
    }

    //oriented

    #[derive(Debug)]
    pub struct Solution{
        pub id_a : usize,
        pub id_b : usize,
        pub orientation : Orientation,
        pub overlap_a : usize,
        pub overlap_b : usize,
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
        pub text : Vec<u8>,
//        pub id2string_vec : Vec<&'a [u8]>,
        pub id2name_vec : Vec<String>,
        pub id2index_bdmap : BidirMap<usize, usize>,
        pub num_ids : usize,
    }

    impl Maps{
        pub fn get_string(&self, id : usize) -> &[u8]{
            assert!(id < self.num_ids);
            &self.text[*self.id2index_bdmap.get_by_first(&id).expect("GAH")..self.get_end_index(id)]
        }

        pub fn get_length(&self, id : usize) -> usize{
            assert!(id < self.num_ids);
            self.get_end_index(id) - self.id2index_bdmap.get_by_first(&id).expect("WOO")
        }

        fn get_end_index(&self, id : usize) -> usize{
            assert!(id < self.num_ids);
            if id == self.num_ids-1{
                self.text.len() - 1 //$s in front. one # at the end
            }else{
                self.id2index_bdmap.get_by_first(&(id + 1)).expect("WAHEY") - 1
            }
        }

        pub fn find_id_for_index_within(&self, index : usize) -> usize{
            let mut best_id = 1;
            for &(id, ind) in self.id2index_bdmap.iter(){
                if index <= ind && ind > best_id{
                    best_id = id;
                }
            }
            best_id
        }

        pub fn print_text_debug(&self){
            println!("{}", String::from_utf8_lossy(&self.text));
        }
    }

    #[derive(Debug)]
    pub struct Config{
        //TODO benchmark argument

        //required
        pub input : String,
        pub output : String,
        pub err_rate : f32,
        pub thresh : i32,
        pub worker_threads: usize,

        //optional
        pub sorted : bool,
        pub reversals : bool,
        pub inclusions : bool,
        pub edit_distance : bool,
        pub verbose : bool,
        pub time: bool,
        pub print: bool,
    }
}

