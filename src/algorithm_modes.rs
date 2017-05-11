


pub mod kucherov{

    use std::cmp::max;

    const S_PARAM : i32 = 3;

    struct Filter{
        blocks : i32,
    }

    impl IntoIterator for Filter {
        type Item = i32;
        type IntoIter = FilterIterator;

        fn into_iter(self) -> Self::IntoIter {
            FilterIterator { blocks: self.blocks, next_index: 0 }
        }
    }

    pub struct FilterIterator{
        blocks : i32,
        next_index : i32,
    }

    impl Iterator for FilterIterator {
        type Item = i32;
        fn next(&mut self) -> Option<i32> {
            if self.next_index == self.blocks{
                None
            } else{
                self.next_index += 1;
                Some(if self.next_index == self.blocks {self.next_index-2} else {self.next_index-1})
            }
        }
    }

    pub fn block_lengths_and_filters(patt_len : i32, err_rate : f32, thresh : i32)
                                    -> (Vec<i32>, Vec<FilterIterator>){
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

        println!("{:?} blocks!", block_lengths.len());

        let mut filters : Vec<FilterIterator> = Vec::new();
        for i in (S_PARAM-1..(block_lengths.len()) as i32).rev() {
            filters.push(FilterIterator{blocks:(i+1) as i32, next_index:0});
        }

        (block_lengths, filters)
    }

    pub fn candidate_condition(
                p_i_start : i32,
                p_i_next  : i32,
                thresh    : i32,
                block_id  : i32,
                errors    : i32
                ) -> bool{
        let c1 = p_i_next >= thresh;
        let c2 = block_id > 0;
        let c3 = block_id >= S_PARAM-1
            && errors <= (block_id - S_PARAM + 1);
        c1 && c2 && c3
    }
}