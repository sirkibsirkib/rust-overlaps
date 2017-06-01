
/*
Kucherov's suffix filters and partition scheme for the ASPOP exact algorithm.
This module is loaded so that it's functions can be used in the candidate generation step.

Other modules can be swapped in as long as they have the same functions.
Ensure that your scheme is CORRECT (doesn't miss solutions)
TODO write tests to check
*/



pub mod kucherov{

    use std::cmp::{min, max};

    pub const S_PARAM : i32 = 2;

    #[inline]
    pub fn get_guaranteed_extra_blocks() -> i32 {
        S_PARAM
    }

    //[~~~~suff blocks~~~|~~~~~blinc blocks~~~~]
    // patt block index is which block we are in relative to PATT
    //  ie: [5,4,3,2,1,0] for all suffs
    #[inline]
    pub fn filter_func(completed_blocks : i32, patt_blocks : i32, blind_blocks : i32) -> i32{
        min(
            completed_blocks,
            patt_blocks - S_PARAM,
        )
    }

    #[inline]
    pub fn get_block_lengths(patt_len : i32, err_rate : f32, thresh : i32) -> Vec<i32>{
        let mut ls : Vec<i32> = Vec::new();
        for l in thresh..patt_len+1{
            let f_len = l as f32;
            if (err_rate*(f_len-1.0)).ceil() < (err_rate*f_len).ceil() {
                ls.push(l);
            }
        }
        ls.push(patt_len+1);
        let k : i32 = (err_rate*(ls[0] as f32) + (S_PARAM as f32) - 1.0).ceil() as i32;
        let a : i32 = (((ls[0]-1) as f32)/(k as f32)).ceil() as i32;
        let b : i32 = ls[0]-thresh;
        let l : i32 = max(a, b);
        let p : i32 =  ((ls[0]-1-l) as f32 / ((k-1) as f32)).floor() as i32;
        let first_half_len : i32 = p*(k-1)+l;
        let longer_blocks_in_first_half : i32 = ls[0]-1-first_half_len;

        let mut block_lengths = Vec::new();
        for _ in 0..(k-1-longer_blocks_in_first_half){
            block_lengths.push(p);
        }
        for _ in 0..longer_blocks_in_first_half{
            block_lengths.push(p+1);
        }
        block_lengths.push(l);
        for i in 0..ls.len()-1 {
            block_lengths.push(ls[i+1] - ls[i]);
        }
        block_lengths
    }

    /*
    During the recursive index search of the candidate generation step, this function will be called
    to check whether it is necessary to generate candidates at this point.
    NOTE: Be more lenient toward correctness
    NOTE: Be more strict toward a faster candidate verification phase
    */
    #[inline]
    pub fn candidate_condition(
                generous_match_len : i32,
                completed_blocks : i32,
                thresh : i32,
                errors : i32
                ) -> bool{
        let c1 = generous_match_len >= thresh;
        let c2 = completed_blocks > 0;
        let c3 = completed_blocks >= S_PARAM - 1
            && errors <= (completed_blocks - S_PARAM + 1);
        c1 && c2 && c3
    }
}


pub mod valimaki2{
    use std::cmp::{min, max};

    #[inline]
    pub fn get_guaranteed_extra_blocks() -> i32 {
        1
    }

    #[inline]
    pub fn filter_func(completed_blocks : i32, patt_blocks : i32, blind_blocks : i32) -> i32{
        if blind_blocks == 0 {
            completed_blocks + 1
        } else {
            completed_blocks
        }
    }

    #[inline]
    pub fn get_block_lengths(patt_len : i32, err_rate : f32, thresh : i32) -> Vec<i32>{
        let mut ps : Vec<i32> = Vec::new();
        assert!(thresh <= patt_len);
        for l in thresh..patt_len+1{
            let one_p = (
                (l as f32)
                    /
                ((err_rate * (l as f32)).ceil() + 1.0)
            ).ceil() as i32;
            ps.push(one_p);
        }
        println!("{:?}", &ps);
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
    pub fn candidate_condition(
            generous_match_len : i32,
            completed_blocks : i32,
            thresh : i32,
            errors : i32
            ) -> bool{
        let c1 = generous_match_len >= thresh;
        let c2 = completed_blocks > 0;
        c1 && c2
    }
}




///////////////// KUCHEROV ///////////////
//#[inline]
//#[cfg(feature="kucherov")]
//pub fn get_guaranteed_extra_blocks() -> i32 {
//    kucherov::get_guaranteed_extra_blocks()
//}
//
//
//#[inline]
//#[cfg(feature="kucherov")]
//pub fn get_fewest_suff_blocks() -> i32{
//    kucherov::fewest_suff_blocks()
//}
//
//#[inline]
//#[cfg(feature="kucherov")]
//pub fn filter_func(completed_blocks : i32, patt_blocks : i32, blind_blocks : i32) -> i32{
//    kucherov::filter_func(completed_blocks, patt_blocks, blind_blocks)
//}
//
//#[inline]
//#[cfg(feature="kucherov")]
//pub fn get_block_lengths(patt_len : i32, err_rate : f32, thresh : i32) -> Vec<i32>{
//    kucherov::get_block_lengths(patt_len, err_rate, thresh)
//}
//
//#[inline]
//#[cfg(feature="kucherov")]
//pub fn candidate_condition(generous_overlap_len : i32, completed_blocks : i32, thresh : i32, errors : i32 ) -> bool{
//    kucherov::candidate_condition(generous_overlap_len, completed_blocks, thresh, errors)
//}
//
///////////////// VALIMAKI ///////////////
//
//#[inline]
//#[cfg(feature="valimaki2")]
//pub fn get_guaranteed_extra_blocks() -> i32 {
//    valimaki2::get_guaranteed_extra_blocks()
//}
//
//
//#[inline]
//#[cfg(feature="valimaki2")]
//pub fn get_fewest_suff_blocks() -> i32{
//    valimaki2::fewest_suff_blocks()
//}
//
//#[inline]
//#[cfg(feature="valimaki2")]
//pub fn filter_func(completed_blocks : i32, patt_blocks : i32, blind_blocks : i32) -> i32{
//    valimaki2::filter_func(completed_blocks, patt_blocks, blind_blocks)
//}
//
//
//#[inline]
//#[cfg(feature="valimaki2")]
//pub fn get_block_lengths(patt_len : i32, err_rate : f32, thresh : i32) -> Vec<i32>{
//    valimaki2::get_block_lengths(patt_len, err_rate, thresh)
//}
//
//#[inline]
//#[cfg(feature="valimaki2")]
//pub fn candidate_condition(generous_overlap_len : i32, completed_blocks : i32, thresh : i32, errors : i32 ) -> bool{
//    valimaki2::candidate_condition(generous_overlap_len, completed_blocks, thresh, errors)
//}




//////////////////////////////////////////////////////////////////
// used only by --tests
#[inline]
pub fn get_guaranteed_extra_blocks() -> i32 {
    S_PARAM
}

pub fn fewest_suff_blocks() -> i32{
    S_PARAM
}

//[~~~~suff blocks~~~|~~~~~blind blocks~~~~]
// patt block index is which block we are in relative to PATT
//  ie: [5,4,3,2,1,0] for all suffs
#[inline]
pub fn filter_func(completed_blocks : i32, patt_blocks : i32, blind_blocks : i32) -> i32{
    min(
        completed_blocks,
        patt_blocks - S_PARAM,
    )
}

#[inline]
pub fn get_block_lengths(patt_len : i32, err_rate : f32, thresh : i32) -> Vec<i32>{
    let mut ls : Vec<i32> = Vec::new();
    for l in thresh..patt_len+1{
        let f_len = l as f32;
        if (err_rate*(f_len-1.0)).ceil()
            < (err_rate*f_len).ceil() {
            ls.push(l);
        }
    }
    ls.push(patt_len+1);
    let k = (err_rate*(ls[0] as f32)).ceil() as i32 + S_PARAM - 1;
    let big_l : i32 = max(
        (((ls[0]-1) as f32)/(k as f32)).ceil() as i32,
        ls[0] - thresh,
    );
    let p : i32 =  ((ls[0]-1-big_l) as f32 / ((k-1) as f32)).floor() as i32;
    let first_half_len : i32 = p*(k-1)+big_l;
    let longer_blocks_in_first_half : i32 = ls[0]-1-first_half_len;

    let mut block_lengths = Vec::new();
    for _ in 0..(k-1-longer_blocks_in_first_half){
        //shorter PRIOR blocks
        block_lengths.push(p);
    }
    for _ in 0..longer_blocks_in_first_half{
        //longer prior blocks
        block_lengths.push(p+1);
    }
    //L block
    block_lengths.push(big_l);
    for i in 0..ls.len()-1 {
        //ANTERIOR blocks
        block_lengths.push(ls[i+1] - ls[i]);
    }
    assert_eq!(block_lengths[0..(k as usize)].iter().sum::<i32>(), ls[0] - 1);
    assert!(block_lengths[(k-1) as usize] >= ls[0] - thresh);
    block_lengths
}

//pub fn get_block_lengths(patt_len : i32, err_rate : f32, thresh : i32) -> Vec<i32>{
//    let mut ls : Vec<i32> = Vec::new();
//    for l in thresh..patt_len+1{
//        let f_len = l as f32;
//        if (err_rate*(f_len-1.0)).ceil() < (err_rate*f_len).ceil() {
//            ls.push(l);
//        }
//    }
//    ls.push(patt_len+1);
//    let k : i32 = (err_rate*(ls[0] as f32) + (S_PARAM as f32) - 1.0).ceil() as i32;
//    let a : i32 = (((ls[0]-1) as f32)/(k as f32)).ceil() as i32;
//    let b : i32 = ls[0]-thresh;
//    let l : i32 = max(a, b);
//    let p : i32 =  ((ls[0]-1-l) as f32 / ((k-1) as f32)).floor() as i32;
//    let first_half_len : i32 = p*(k-1)+l;
//    let longer_blocks_in_first_half : i32 = ls[0]-1-first_half_len;
//
//    let mut block_lengths = Vec::new();
//    for _ in 0..(k-1-longer_blocks_in_first_half){
//        //shorter PRIOR blocks
//        block_lengths.push(p);
//    }
//    for _ in 0..longer_blocks_in_first_half{
//        //longer prior blocks
//        block_lengths.push(p+1);
//    }
//    //L block
//    block_lengths.push(l);
//    for i in 0..ls.len()-1 {
//        //ANTERIOR blocks
//        block_lengths.push(ls[i+1] - ls[i]);
//    }
//    block_lengths
//}

/*
During the recursive index search of the candidate generation step, this function will be called
to check whether it is necessary to generate candidates at this point.
NOTE: Be more lenient toward correctness
NOTE: Be more strict toward a faster candidate verification phase
*/
#[inline]
pub fn candidate_condition(
    generous_overlap_len : i32,
    completed_blocks : i32,
    thresh : i32,
    errors : i32
) -> bool{
    let c1 = generous_overlap_len >= thresh;
    let c2 = completed_blocks > 0;
    let c3 = completed_blocks >= S_PARAM - 1
        &&
        errors <= (completed_blocks - S_PARAM + 1);
    c1 && c2 && c3
}