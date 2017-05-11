pub mod solutions{

    pub enum Orientation{
        Normal,
        Inverted,
    }

    //should already be correctly oriented
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
}