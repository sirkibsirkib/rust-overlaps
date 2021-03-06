use std::cmp::{min};
use modes::IsMode;
use std::fmt;

#[derive(Debug)]
pub struct ValimakiMode;

impl ValimakiMode {
    pub fn new() -> Self{
        ValimakiMode
    }
}

/*
According to the Valimaki algorithm in the 2012 paper
*/
impl fmt::Display for ValimakiMode {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Valimaki")
    }
}

#[allow(unused_variables)]
impl IsMode for ValimakiMode {

    #[inline]
    fn get_guaranteed_extra_blocks(&self) -> i32 {
        1
    }

    fn get_fewest_suff_blocks(&self) -> i32{
        1
    }

    #[inline]
    fn filter_func(&self, completed_blocks : i32, patt_blocks : i32, blind_blocks : i32) -> i32{
        if blind_blocks == 0 {
            completed_blocks + 1 + 1
        } else {
            completed_blocks + 1
        }
    }

    #[inline]
    fn get_block_lengths(&self, patt_len : i32, err_rate : f32, thresh : i32) -> Vec<i32>{
        let mut ps : Vec<i32> = Vec::new();
        if patt_len < thresh{
            ps.push(patt_len);
            return ps;
        }
        for l in thresh..patt_len+1{
            let one_p = (
                (l as f32)
                    /
                    ((err_rate * (l as f32)).ceil() + 1.0)
            ).ceil() as i32;
            ps.push(one_p);
        }
        let p = *ps.iter().min().unwrap();
        let mut remain = patt_len;
        let mut block_lengths : Vec<i32> = Vec::new();
        while remain > 0{
            let next = min(remain, p);
            block_lengths.push(next);
            remain -= next;
        }
        block_lengths
    }

    #[inline]
    fn candidate_condition(&self,
        generous_overlap_len : i32,
        completed_blocks : i32,
        thresh : i32,
        errors : i32
    ) -> bool{
        let c1 = generous_overlap_len >= thresh;
        let c2 = completed_blocks > 0;
        c1 && c2
    }
}