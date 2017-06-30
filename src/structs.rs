use useful;

/*
Structs and methods for solutions and candidate solutions to be used throughout the program
*/
pub mod solutions{
    use std::hash::{Hash, Hasher};
    use std::cmp::Ordering;
    use std::cmp::max;
    use std::mem::swap;
    use super::useful::{companion_id, Orientation};

    //NOT oriented
    #[derive(Hash,PartialEq, Eq, Debug, Clone)]
    pub struct Candidate{
        pub id_b : usize,
        pub overlap_a : usize,
        pub overlap_b : usize,
        pub overhang_left_a : i32,

        //DEBUG
//        pub debug_str : String,
    }

    impl Candidate{

        #[inline]
        pub fn a1(&self) -> usize {
            max(0, self.overhang_left_a) as usize
        }

        #[inline]
        pub fn b1(&self) -> usize {
            max(0, -self.overhang_left_a) as usize
        }

        #[inline]
        pub fn a2(&self) -> usize {
            self.overlap_a
        }

        #[inline]
        pub fn b2(&self) -> usize {
            self.overlap_b
        }

        //TODO tidy up these silly a2() etc. calls later
        #[inline]
        pub fn a3(&self, a_len : usize) -> usize {
            assert!(a_len >= self.a1() + self.overlap_a);
            a_len - self.a1() - self.a2()
        }

        #[inline]
        pub fn b3(&self, b_len : usize) -> usize {
            assert!(b_len >= self.b1() + self.overlap_b);
            b_len - self.b1() - self.b2()
        }
    }

    //oriented
    #[derive(Debug,Clone)]
    pub struct Solution{
        pub id_a : usize,
        pub id_b : usize,
        pub orientation : Orientation,
        pub overhang_left_a : i32,
        pub overhang_right_b : i32,
        pub overlap_a : usize,
        pub overlap_b : usize,
        pub errors : u32,
    }

    impl Solution{
        pub fn v_flip(&mut self){
            self.overhang_left_a *= -1;
            self.overhang_right_b *= -1;
            swap(&mut self.id_a, &mut self.id_b);
            swap(&mut self.overlap_a, &mut self.overlap_b);
        }

        pub fn h_flip(&mut self, reversals : bool){
            self.id_a = companion_id(self.id_a, reversals);
            self.id_b = companion_id(self.id_b, reversals);
            self.mirror_horizontally();
        }

        //strictly reverses orientation to compensate for index being backwards
        pub fn mirror_horizontally(&mut self){
            swap(&mut self.overhang_left_a, &mut self.overhang_right_b);
            self.overhang_left_a *= -1;
            self.overhang_right_b *= -1;
        }
    }

    impl Ord for Solution {
        fn cmp(&self, other: &Self) -> Ordering {
            (self.id_a, self.id_b, &self.orientation, self.overhang_left_a, self.overhang_right_b, self.overlap_a, self.overlap_b)
                .cmp(&(other.id_a, other.id_b, &other.orientation, other.overhang_left_a, other.overhang_right_b, other.overlap_a, other.overlap_b))
        }
    }

    impl PartialOrd for Solution {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            Some(self.cmp(other))
        }
    }

    impl PartialEq for Solution {
        fn eq(&self, other: &Self) -> bool {
            (self.id_a, self.id_b, &self.orientation, self.overhang_left_a, self.overhang_right_b, self.overlap_a, self.overlap_b)
                == (other.id_a, other.id_b, &other.orientation, other.overhang_left_a, other.overhang_right_b, other.overlap_a, other.overlap_b)
        }
    }

    impl Eq for Solution { }

    impl Hash for Solution {
        fn hash<H: Hasher>(&self, state: &mut H) {
            self.id_a.hash(state);
            self.id_b.hash(state);
            self.orientation.hash(state);
            self.overlap_a.hash(state);
            self.overlap_b.hash(state);
            self.overhang_left_a.hash(state);
            self.overhang_right_b.hash(state);
        }
    }
}


/*
Some structs and convenience functions for storing the data needed throughout the run.
The config struct stores the user's input parameters and is checked frequently but never changed.
The maps struct is built once, stores text, and various mappings between different string
representations and is queried throughout the run, also never changing after being populated.
*/
pub mod run_config{
    extern crate bidir_map;
    use bidir_map::BidirMap;

    #[derive(Debug)]
    pub struct Maps{
        pub text : Vec<u8>,
        pub id2name_vec : Vec<String>,
        pub id2index_bdmap : BidirMap<usize, usize>,
        pub indexes : Vec<usize>,
    }

    impl Maps{

        pub fn num_ids(&self) -> usize {
            self.id2index_bdmap.len()
        }

        pub fn get_string(&self, id : usize) -> &[u8]{
            assert!(id < self.num_ids());
            &self.text[*self.id2index_bdmap.get_by_first(&id).unwrap()..self.get_end_index(id)]
        }

        pub fn get_length(&self, id : usize) -> usize{
            assert!(id < self.num_ids());
            self.get_end_index(id) - self.id2index_bdmap.get_by_first(&id).unwrap()
        }

        fn get_end_index(&self, id : usize) -> usize{
            assert!(id < self.num_ids());
            if id == self.num_ids()-1{
                self.text.len() - 1 //$s in front. one # at the end
            }else{
                self.id2index_bdmap.get_by_first(&(id + 1)).unwrap() - 1
            }
        }

        //returns (id, index)
        #[inline]
        pub fn find_occurrence_containing(&self, index : usize) -> (usize, usize){
            match self.indexes.binary_search(&index){
                Ok(found_id) => (found_id, index),
                Err(insert_id) => (insert_id-1, self.index_for(insert_id-1)),
            }
        }

        pub fn get_name_for(&self, id : usize) -> &str {
            self.id2name_vec.get(id).expect("get name")
        }

        #[inline]
        pub fn id_for(&self, id : usize) -> usize{
            *(self.id2index_bdmap.get_by_second(&id)
                .expect(&format!("no index for ID {}. input has IDs from 0 --> {}",
                                id, self.num_ids())))
        }

        #[inline]
        pub fn index_for(&self, index : usize) -> usize{
            *(self.id2index_bdmap.get_by_first(&index)
                .expect(&format!("no id at index {}", index)))
        }
    }

    pub static N_ALPH : &'static [u8] = b"ACGNT";
    pub static ALPH : &'static [u8] = b"ACGT";


    #[derive(Debug)]
    pub struct Config{

        //required
        pub input : String,
        pub output : String,
        pub err_rate : f32,
        pub thresh : i32,

        //optional
        pub format_line: bool,
        pub greedy_output: bool,
        pub reversals : bool,
        pub inclusions : bool,
        pub edit_distance : bool,
        pub verbosity: u8,
        pub print: bool,
        pub n_alphabet: bool,
        pub track_progress: bool,
        pub worker_threads: usize,
    }

    impl Config{
        pub fn alphabet(&self) -> &[u8]{
            if self.n_alphabet {
                &N_ALPH
            } else {
                &ALPH
            }
        }
    }
}
