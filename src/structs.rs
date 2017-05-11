pub mod solutions{
    pub enum Orientation{
        Normal,
        Inverted,
    }

    //should already be correctly oriented
    #[derive(Hash)]
    pub struct Candidate{
        id_b : i32,
        orientation : Orientation,
        overlap_a : i32,
        overlap_b : i32,
        overhang_left_a : i32,
        overhang_right_b : i32,
    }

    #[derive(Hash)]
    pub struct Solution{
        id_b : i32,
        orientation : Orientation,
        overlap_a : i32,
        overlap_b : i32,
        overhang_left_a : i32,
        overhang_right_b : i32,
        errors : u8,
        cigar : String,
    }

    impl Hash for Struct {
        fn hash<H: Hasher>(&self, state: &mut H) {
            self.errors.hash(state);
            self.candidate.hash(state);
            // cigar string not involved in the hash
        }
    }
}

pub mod run_config{
    #[derive(Debug)]
    pub struct Maps{
        id2name : HashMap<i32, String>,
        id2str_in_s : HashMap<i32, Vec<u8>>,
        bdmap_index_id : BidirMap<i32, i32>,
    }
    impl Maps{
        pub fn num_strings(&self) -> usize{
            self.id2str_in_s.len()
        }
    }

    #[derive(Debug)]
    pub struct Config{
        //required
        input : String,
        output : String,
        err_rate : f32,
        thresh : i32,
        worker_threads: i32,

        //optional
        reversals : bool,
        inclusions : bool,
        edit_distance : bool,
        verbose : bool,
    }
}