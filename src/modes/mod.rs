mod kucherov;
mod valimaki2;

/////////////// KUCHEROV ///////////////
#[inline]
#[cfg(feature="kucherov")]
pub fn get_guaranteed_extra_blocks() -> i32 {
    kucherov::get_guaranteed_extra_blocks()
}

#[inline]
#[cfg(feature="kucherov")]
pub fn filter_func(completed_blocks : i32, patt_blocks : i32, blind_blocks : i32) -> i32{
    kucherov::filter_func(completed_blocks, patt_blocks, blind_blocks)
}


#[inline]
#[cfg(feature="kucherov")]
pub fn get_block_lengths(patt_len : i32, err_rate : f32, thresh : i32) -> Vec<i32>{
    kucherov::get_block_lengths(patt_len, err_rate, thresh)
}

#[inline]
#[cfg(feature="kucherov")]
pub fn candidate_condition(generous_match_len : i32, completed_blocks : i32, thresh : i32, errors : i32 ) -> bool{
    kucherov::candidate_condition(generous_match_len, completed_blocks, thresh, errors)
}

/////////////// VALIMAKI ///////////////

#[inline]
#[cfg(feature="valimaki2")]
pub fn get_guaranteed_extra_blocks() -> i32 {
    valimaki2::get_guaranteed_extra_blocks()
}

#[inline]
#[cfg(feature="valimaki2")]
pub fn filter_func(completed_blocks : i32, patt_blocks : i32, blind_blocks : i32) -> i32{
    valimaki2::filter_func(completed_blocks, patt_blocks, blind_blocks)
}


#[inline]
#[cfg(feature="valimaki2")]
pub fn get_block_lengths(patt_len : i32, err_rate : f32, thresh : i32) -> Vec<i32>{
    valimaki2::get_block_lengths(patt_len, err_rate, thresh)
}

#[inline]
#[cfg(feature="valimaki2")]
pub fn candidate_condition(generous_match_len : i32, completed_blocks : i32, thresh : i32, errors : i32 ) -> bool{
    valimaki2::candidate_condition(generous_match_len, completed_blocks, thresh, errors)
}