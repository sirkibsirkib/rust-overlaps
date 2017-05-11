pub mod solutions{

    #[derive(Hash)]
    pub struct RawCandidate{
        id_b : i32,
        overlap_a : i32,
        overlap_b : i32,
        overhang_right_b : i32,
    }


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

    pub struct Solution<'a>{
        candidate : &'a Candidate,
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