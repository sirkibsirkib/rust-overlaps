
use structs::solutions::*;
use structs::run_config::*;


fn step2(text : &str, config : &Config, maps : &Maps){
    //TODO thread_pool_2
    for p_id in maps.id2str_in_s.keys(){
        let candidates : HashSet<Candidate> = read_candidates(p_id);
        for rc in candidates{

        }
    }
}


fn read_candidates(p_id : i32) -> HashSet<Candidate>{
    let result : HashSet<Candidate> = HashSet::new();
    let filepath = TEMP_DIR.to_path_buf().push(p_id.to_string());
    let mut rdr = csv::Reader::from_file(Path::new(filepath));
    for row in rdr.decode(){
        let (id_other, overlap_a, overlap_b, overhang_right_b): (i32, i32, u32, i32) = row.unwrap();
        result.
    }
    //flush writer? idk
}

fn test(){

}


fn verify_all(candidates : &HashSet<Candidate>) -> HashSet<Solution>{
    let solution_set : HashSet<Solution> = HashSet:new();
    for c in candidates{
        //TODO
        match verify()
    }
}

// incoming cands are already oriented
fn verify(id_a : i32, c : &Candidate, config : &Config, maps : &Maps) -> Option(Solution){
    let a_len = maps.id2str_in_s.get(id_a).len();
    let b_len = maps.id2str_in_s.get(c.id_b).len();


    let a_start = max(0, c.overhang_left_a);
    let a_end = a_start + c.overlap_a;
    let a_end2 = min(a_len, a_len-c.overhang_right_b);
    assert_eq!(a_end, a_end2);          //DEBUG! can remove a_end2 if this checks out

    let b_start = max(0, -c.overhang_left_a);
    let b_end = b_start + c.overlap_b;
    let b_end2 = min(b_len, b_len+c.overhang_right_b);
    assert_eq!(b_end, b_end2);          //DEBUG! can remove b_end2 if this checks out

    let a_part : &str = maps.id2str_in_s.get(id_a)  [a_start..a_end];
    let b_part : &str = maps.id2str_in_s.get(c.id_b)[b_start..b_end];

    let k : i32 = if config.edit_distance{
        hamming(a_part, b_part).expect("Can't apply hamming to mismatching lengths!")
    }else{
        levenshtein(a_part, b_part)
    };
    let k_limit : f32 = config.err_rate*(max(c.overlap_a, c.overlap_b) as f32);
    if k <= k_limit{
        Some(Solution{candidate:c, errors:k, cigar:""})
    }else{
        None()
    }
}