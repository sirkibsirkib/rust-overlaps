
use structs::solutions::*;
use structs::run_config::*;

fn step2(text : &str, config : &Config, maps : &Maps) -> Vec<HashSet<Candidate>>{
    let p_ids = maps.id2str_in_s.keys();
    let mut all_solutions : Vec<HashSet<Candidate>> = HashSet::with_capacity(p_ids.len());
    for p_id in p_ids{
        //TODO thread_pool_2
        let candidates : HashSet<Candidate> = read_candidates(p_id);
        all_solutions[p_id] = (verify_all(p_id, candidates);
    }
    all_solutions
}

// figure out boxes, seriously;
fn read_candidates(p_id : usize) -> Box<HashSet<Candidate>>{
    let result : Box<HashSet<Candidate>> = Box::new(HashSet::new());
    let filepath = TEMP_DIR.to_path_buf().push(p_id.to_string());
    let mut rdr = csv::Reader::from_file(Path::new(filepath));
    for row in rdr.decode() {
        let (id_other, overlap_a, overlap_b, overhang_right_b): (usize, i32, u32, i32) = row.unwrap();
        let candidate = Candidate {
            id_b: id_other,
            overlap_a: overlap_a,
            overlap_b: overlap_b,
            overhang_right_b: overhang_right_b,
        };
        result.insert(candidate)
    }
    result
}


fn verify_all(p_id : usize, candidates : &HashSet<Candidate>) -> HashSet<Solution>{
    let solution_set : HashSet<Solution> = HashSet:new();
    for c in candidates{
        if let Some(solution) = solution_set.insert(solution){
            solution_set.insert(solution);
        }
    }
    solution_set
}

// incoming cands are already oriented
fn verify(id_a : i32, c : &Candidate, config : &Config, maps : &Maps) -> Option(Solution){
    let a_len = maps.id2str_in_s.get(id_a).len();
    let b_len = maps.id2str_in_s.get(c.id_b).len();

    let a_start = max(0, c.overhang_left_a);
    let a_end = a_start + c.overlap_a;
    let a_end2 = min(a_len, a_len-c.overhang_right_b);
    //TODO return None here if you detect you're out of bounds of one of the strings
    assert_eq!(a_end, a_end2);          //DEBUG! can remove a_end2 if this checks out

    let b_start = max(0, -c.overhang_left_a);
    let b_end = b_start + c.overlap_b;
    let b_end2 = min(b_len, b_len+c.overhang_right_b);
    //TODO return None here if you detect you're out of bounds of one of the strings
    assert_eq!(b_end, b_end2);          //DEBUG! can remove b_end2 if this checks out

    let a_part : &str = maps.id2str_in_s.get(id_a)  [a_start..a_end];
    let b_part : &str = maps.id2str_in_s.get(c.id_b)[b_start..b_end];

    let errors : i32 = if config.edit_distance{
        hamming(a_part, b_part).expect("Can't apply hamming to mismatching lengths!")
    }else{
        levenshtein(a_part, b_part)
    };
    let k_limit : f32 = config.err_rate*(max(c.overlap_a, c.overlap_b) as f32);
    if errors <= k_limit{
        Some(solution_from_candidate(c, id_a, cigar, errors, maps))
    }else{
        None()
    }
}

fn solution_from_candidate(c : &Candidate, id_a : usize, cigar : String, errors : i32, maps : &Maps) -> Solution {
    let a_len = maps.get(id_a).expect("WTF").len();
    let b_len = maps.get(c.other_id).expect("WTF").len();
    let id_b = c.id_b;

    let overlap_a = c.overlap_a;
    let overlap_b = c.overlap_b;
    let (overhang_a, overhang_b) = if c.overhang_right_b == 0 {
        (
            -(b_len - overlap_b),
            -(a_len - overlap_a),
        )
    } else {
        (
            -(b_len - c.overhang_right_b - overlap_b),
            c.overhang_right_b,
        )
    };
    let (overhang_left_a, overhang_right_b) = if c.overhang_right_b == 0 {
        (
            -(b_len - overlap_b),
            -(a_len - overlap_a),
        )
    } else {
        (
            -(b_len - c.overhang_right_b - overlap_b),
            c.overhang_right_b,
        )
    };

    //START TRANSFORM

    //ensure smaller string comes first
    if a_id > b_id {
        //
        overhang_left_a *= -1;
        overhang_right_b = -1;
        (a_len, b_len) = (b_len, a_len);
        (id_a, id_b) = (id_b, id_a);
        (overlap_a, overlap_b) = (overlap_b, overlap_a);
        cigar = cigar.vflip();  // I->D, D->I
    }

    //first string is a reverse
    if maps.reverses && id_a % 2 == 1 {
        (overhang_a, overhang_b) = (-overhang_b, -overhang_a);
        a_id = companion_id(a_id);
        b_id = companion_id(b_id);
        cigar.h_flip(); // XYZ -> ZYX
    }

    let orientation = if !config.reverses || id_b%2==0{
        structs::Orientation::Normal
    }else{
        structs::Orientation::Reversed
    };
    let a_name = maps.id2name.get(id_a);
    let b_name = maps.id2name.get(id_b);
    structs::Solution{
        id_a : id_a,
        id_b : id_b,
        orientation : orientation,
        overlap_a : overlap_a,
        overlap_b : overlap_b,
        overhang_left_a : overhang_left_a,
        overhang_right_b : overhang_right_b,
        errors : errors,
        cigar : cigar,
    }
}

fn companion_id(id : usize) -> usize{
    if id%2==0 {id+1} else {id-1}
}