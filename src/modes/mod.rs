pub mod kucherov;
pub mod valimaki2;

pub type Mode = Box<IsMode>;

pub trait IsMode: Sync {
    fn get_guaranteed_extra_blocks(&self) -> i32;
    fn get_fewest_suff_blocks(&self) -> i32;
    fn filter_func(&self, completed_blocks : i32, patt_blocks : i32, blind_blocks : i32) -> i32;
    fn get_block_lengths(&self, patt_len : i32, err_rate : f32, thresh : i32) -> Vec<i32>;
    fn candidate_condition(&self,generous_overlap_len : i32, completed_blocks : i32, thresh : i32, errors : i32 ) -> bool;
}

