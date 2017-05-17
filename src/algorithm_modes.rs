pub fn dank() {println!("LOL");}

#[allow(dead_code)]
#[allow(unused_variables)]
pub mod kucherov{

    use std::cmp::{min, max};

    const S_PARAM : i32 = 2;

    //[~~~~suff blocks~~~|~~~~~blinc blocks~~~~]
    // patt block index is which block we are in relative to PATT
    //  ie: [5,4,3,2,1,0] for all suffs
    pub fn filter_func(completed_blocks : i32, patt_blocks : i32) -> i32{
        let x = min(
            completed_blocks,
            patt_blocks - S_PARAM,
        );
//        println!("({}, {}) ===> ({})", completed_blocks, patt_blocks, x);
        x
    }

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