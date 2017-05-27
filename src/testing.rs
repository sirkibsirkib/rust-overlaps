/*
Some simple tests to check if the system and all its various modes behave correctly and find
the right solutions, using some predetermined inputs in th text_input directory.
*/
#[cfg(test)]
mod tests {
    use super::super::*;
    use useful::Orientation;
    use useful::Orientation::*;
    use std::io::{BufReader, BufRead};

    #[derive (Eq, PartialEq, Hash, Debug)]
    struct GoodSolution{
        a_nm: String,
        b_nm: String,
        or: Orientation,
        oha: i32,
        ohb: i32,
        ola: usize,
        olb: usize,
        err: u32,
    }

    #[test]
    fn basic_mapping() {
        let config = Config{
            input  :        "./test_input/basic_mapping.fasta".to_owned(),
            output  :       "./test_output/basic_mapping.txt".to_owned(),
            err_rate :      0.03,
            thresh :        8,
            worker_threads: 1,
            reversals :     false,
            inclusions :    false,
            edit_distance : false,
            verbose :       false,
            greedy_output:  false,
            track_progress :false,
            print:          false,
            n_alphabet:     false,
        };
        let maps = prepare::read_and_prepare(&config.input, &config).expect("Couldn't interpret data.");
        assert_eq!(2, maps.num_ids);
        assert_eq!(maps.text.len(), 5+6+1+1+1);
        assert_eq!(maps.id2name_vec[0], "x");
        assert_eq!(maps.id2name_vec[1], "y");
        assert_eq!(maps.get_string(0).len(), 5);
        assert_eq!(maps.get_string(1).len(), 6);
    }

    #[test]
    fn ham() {
        let config = Config{
            input  :        "./test_input/ham.fasta".to_owned(),
            output  :       "./test_output/ham.txt".to_owned(),
            err_rate :      0.02,
            thresh :        4,
            worker_threads: 1,
            reversals :     false,
            inclusions :    false,
            edit_distance : false,
            verbose :       false,
            greedy_output:  false,
            track_progress :false,
            print:          false,
            n_alphabet:     false,
        };
        let maps = prepare::read_and_prepare(&config.input, &config).expect("Couldn't interpret data.");
        solve(&config, &maps);
        let results = read_output(&config.output);
        let mut should_contain : HashSet<GoodSolution> = HashSet::new();
        should_contain.insert(GoodSolution{a_nm:"x".to_owned(), b_nm:"y".to_owned(), or:Normal, oha:5, ohb:10, ola:5, olb:5, err:0});
        panic_if_solutions_missing(results, should_contain);
    }

    #[test]
    fn edit() {
        let config = Config{
            input  :        "./test_input/edit.fasta".to_owned(),
            output  :       "./test_output/edit.txt".to_owned(),
            err_rate :      0.2,
            thresh :        5,
            worker_threads: 1,
            reversals :     false,
            inclusions :    false,
            edit_distance :     true,
            verbose :       false,
            greedy_output:  false,
            track_progress :false,
            print:          false,
            n_alphabet:     false,
        };
        let maps = prepare::read_and_prepare(&config.input, &config).expect("Couldn't interpret data.");
        solve(&config, &maps);
        let results = read_output(&config.output);
        let mut should_contain : HashSet<GoodSolution> = HashSet::new();
        should_contain.insert(GoodSolution{a_nm:"x".to_owned(), b_nm:"y".to_owned(), or:Normal, oha:7, ohb:7, ola:6, olb:7, err:1});
        panic_if_solutions_missing(results, should_contain);
    }

    #[test]
    fn ham_rev() {
        let config = Config{
            input  :        "./test_input/ham_rev.fasta".to_owned(),
            output  :       "./test_output/ham_rev.txt".to_owned(),
            err_rate :      0.02,
            thresh :        5,
            worker_threads: 1,
            reversals :         true,
            inclusions :    false,
            edit_distance : false,
            verbose :       false,
            greedy_output:  false,
            track_progress :false,
            print:          false,
            n_alphabet:     false,
        };
        let maps = prepare::read_and_prepare(&config.input, &config).expect("Couldn't interpret data.");
        solve(&config, &maps);
        let results = read_output(&config.output);
        let mut should_contain : HashSet<GoodSolution> = HashSet::new();
        should_contain.insert(GoodSolution{a_nm:"x".to_owned(), b_nm:"y".to_owned(), or:Reversed, oha:8, ohb:8, ola:8, olb:8, err:0});
        panic_if_solutions_missing(results, should_contain);
    }

    #[test]
    fn ham_incl() {
        let config = Config{
            input  :        "./test_input/ham_incl.fasta".to_owned(),
            output  :       "./test_output/ham_incl.txt".to_owned(),
            err_rate :      0.02,
            thresh :        6,
            worker_threads: 1,
            reversals :     false,
            inclusions :        true,
            edit_distance : false,
            verbose :       false,
            greedy_output:  false,
            track_progress :false,
            print:          false,
            n_alphabet:     false,
        };
        let maps = prepare::read_and_prepare(&config.input, &config).expect("Couldn't interpret data.");
        solve(&config, &maps);
        let results = read_output(&config.output);
        let mut should_contain : HashSet<GoodSolution> = HashSet::new();
        should_contain.insert(GoodSolution{a_nm:"x".to_owned(), b_nm:"y".to_owned(), or:Normal, oha:4, ohb:-4, ola:8, olb:8, err:0});
        panic_if_solutions_missing(results, should_contain);
    }

    #[test]
    fn ham_no_n() {
        let config = Config{
            input  :        "./test_input/ham_no_n.fasta".to_owned(),
            output  :       "./test_output/ham_no_n.txt".to_owned(),
            err_rate :      0.02,
            thresh :        5,
            worker_threads: 1,
            reversals :     false,
            inclusions :    false,
            edit_distance : false,
            verbose :       false,
            greedy_output:  false,
            track_progress :false,
            print:          false,
            n_alphabet:         false,
        };
        let maps = prepare::read_and_prepare(&config.input, &config).expect("Couldn't interpret data.");
        solve(&config, &maps);
        let results = read_output(&config.output);
        let mut should_contain : HashSet<GoodSolution> = HashSet::new();
        should_contain.insert(GoodSolution{a_nm:"x".to_owned(), b_nm:"y".to_owned(), or:Normal, oha:6, ohb:9, ola:6, olb:6, err:0});
        panic_if_solutions_missing(results, should_contain);
    }

    #[test]
    fn ham_rev_incl() {
        let config = Config{
            input  :        "./test_input/ham_rev_incl.fasta".to_owned(),
            output  :       "./test_output/ham_rev_incl.txt".to_owned(),
            err_rate :      0.02,
            thresh :        5,
            worker_threads: 1,
            reversals :         true,
            inclusions :        true,
            edit_distance : false,
            verbose :       false,
            greedy_output:  false,
            track_progress :false,
            print:          false,
            n_alphabet:     false,
        };
        let maps = prepare::read_and_prepare(&config.input, &config).expect("Couldn't interpret data.");
        solve(&config, &maps);
        let results = read_output(&config.output);
        let mut should_contain : HashSet<GoodSolution> = HashSet::new();
        should_contain.insert(GoodSolution{a_nm:"x".to_owned(), b_nm:"y".to_owned(), or:Reversed, oha:5, ohb:-5, ola:5, olb:5, err:0});
        panic_if_solutions_missing(results, should_contain);
    }

    #[test]
    fn edit_rev() {
        let config = Config{
            input  :        "./test_input/edit_rev.fasta".to_owned(),
            output  :       "./test_output/edit_rev.txt".to_owned(),
            err_rate :      0.18,
            thresh :        7,
            worker_threads: 1,
            reversals :         true,
            inclusions :    false,
            edit_distance :     true,
            verbose :       false,
            greedy_output:  false,
            track_progress :false,
            print:          false,
            n_alphabet:     false,
        };
        let maps = prepare::read_and_prepare(&config.input, &config).expect("Couldn't interpret data.");
        solve(&config, &maps);
        let results = read_output(&config.output);
        let mut should_contain : HashSet<GoodSolution> = HashSet::new();
        should_contain.insert(GoodSolution{a_nm:"x".to_owned(), b_nm:"y".to_owned(), or:Reversed, oha:4, ohb:8, ola:7, olb:6, err:1});
        panic_if_solutions_missing(results, should_contain);
    }

    #[test]
    fn edit_incl() {
        let config = Config{
            input  :        "./test_input/edit_incl.fasta".to_owned(),
            output  :       "./test_output/edit_incl.txt".to_owned(),
            err_rate :      0.17,
            thresh :        7,
            worker_threads: 1,
            reversals :     false,
            inclusions :        true,
            edit_distance :     true,
            verbose :       false,
            greedy_output:  false,
            track_progress :false,
            print:          false,
            n_alphabet:     false,
        };
        let maps = prepare::read_and_prepare(&config.input, &config).expect("Couldn't interpret data.");
        solve(&config, &maps);
        let results = read_output(&config.output);
        let mut should_contain : HashSet<GoodSolution> = HashSet::new();
        should_contain.insert(GoodSolution{a_nm:"x".to_owned(), b_nm:"y".to_owned(), or:Normal, oha:4, ohb:-8, ola:7, olb:6, err:1});
        panic_if_solutions_missing(results, should_contain);
    }

    #[test]
    fn edit_rev_incl() {
        let config = Config{
            input  :        "./test_input/edit_rev_incl.fasta".to_owned(),
            output  :       "./test_output/edit_rev_incl.txt".to_owned(),
            err_rate :      0.21,
            thresh :        5,
            worker_threads: 1,
            reversals :         true,
            inclusions :        true,
            edit_distance :     true,
            verbose :       false,
            greedy_output:  false,
            track_progress :false,
            print:          false,
            n_alphabet:     false,
        };
        let maps = prepare::read_and_prepare(&config.input, &config).expect("Couldn't interpret data.");
        solve(&config, &maps);
        let results = read_output(&config.output);
        let mut should_contain : HashSet<GoodSolution> = HashSet::new();
        should_contain.insert(GoodSolution{a_nm:"x".to_owned(), b_nm:"y".to_owned(), or:Reversed, oha:5, ohb:-5, ola:5, olb:6, err:1});
        panic_if_solutions_missing(results, should_contain);
    }

    #[test]
    fn many_errors() {
        let config = Config{
            input  :        "./test_input/many_errors.fasta".to_owned(),
            output  :       "./test_output/many_errors.txt".to_owned(),
            err_rate :      0.4,
            thresh :        8,
            worker_threads: 1,
            reversals :     false,
            inclusions :    false,
            edit_distance : false,
            verbose :       false,
            greedy_output:  false,
            track_progress :false,
            print:          false,
            n_alphabet:     false,

        };
        let maps = prepare::read_and_prepare(&config.input, &config).expect("Couldn't interpret data.");
        solve(&config, &maps);
        let results = read_output(&config.output);
        let mut should_contain : HashSet<GoodSolution> = HashSet::new();
        should_contain.insert(GoodSolution{a_nm:"x".to_owned(), b_nm:"y".to_owned(), or:Normal, oha:0, ohb:0, ola:10, olb:10, err:4});
        panic_if_solutions_missing(results, should_contain);
    }

    #[test]
    fn modified_levenshtein() {
        use verification::modified_levenshtein;

        assert_eq!(modified_levenshtein(b"", b""), 0);
        assert_eq!(modified_levenshtein(b"A", b"A"), 0);
        assert_eq!(modified_levenshtein(b"N", b"N"), 1);
        assert_eq!(modified_levenshtein(b"AA", b"AA"), 0);
        assert_eq!(modified_levenshtein(b"AN", b"AN"), 1);
        assert_eq!(modified_levenshtein(b"AA", b"ATA"), 1);
        assert_eq!(modified_levenshtein(b"ATA", b"AA"), 1);
        assert_eq!(modified_levenshtein(b"AAAAA", b"CAAAC"), 2);
        assert_eq!(modified_levenshtein(b"TTTTA", b"TTTT"), 2);
        assert_eq!(modified_levenshtein(b"G", b""), std::u32::MAX);
        assert_eq!(modified_levenshtein(b"GG", b"G"), std::u32::MAX);
    }

    #[test]
    fn filter_correct() {
//        panic!("VALIMAKI2 not working??try it out");

        use modes::*;
//        use algorithm_modes::valimaki2::*;

        let guaranteed_extra_blocks = get_guaranteed_extra_blocks();
        for patt_len in 5..350 {
            let err_iter = ErrIterator{next:0.5};
            for err_rate in err_iter {
                for thresh in 4..(patt_len as f32 * 0.5)as i32 {
                    let blocks_lengths = get_block_lengths(patt_len, err_rate, thresh);
                    assert!(!(blocks_lengths.contains(&0))); //no-char blocks would mess with the id lookup
                    let mut block_id_lookup = search::get_block_id_lookup(&blocks_lengths);
                    block_id_lookup.reverse();
                    for pref_len in thresh..patt_len+1{
                        let pref_blocks = block_id_lookup[(pref_len-1) as usize] + 1; //#blocks is block index + 1
//                        let effective_pref_blocks = filter_func(pref_blocks-1, blocks_lengths.len() as i32);
                        let max_allowed_err = (pref_len as f32 * err_rate).floor() as i32;

//                        let filter_permits = filter_func(pref_in_blocks, blocks_lengths.len() as i32);
                        if !(pref_blocks >= max_allowed_err + guaranteed_extra_blocks) { // must allow K+1 or more
                            panic!("\nfilter not lenient enough for patt_len {} err_rate {} \
                            thresh {} pref_len {}.\nBlock lens is {:?}. pref in {} blocks, permitted {} errors.\n",
                            patt_len, err_rate, thresh, pref_len, &blocks_lengths, pref_blocks, max_allowed_err);
                        }
                    }
                }
            }
        }
    }

    struct ErrIterator{
        next : f32,
    }

    impl Iterator for ErrIterator {
        type Item = f32;
        fn next(&mut self) -> Option<f32> {
            if self.next < 0.1{
                None
            }else{
                let old_next = self.next;
                self.next *= 0.8;
                Some(old_next)
            }
        }
    }

    fn read_output(filename : &str) -> HashSet<GoodSolution>{
        let mut result : HashSet<GoodSolution> = HashSet::new();
        let f = File::open(filename).unwrap();
        let file = BufReader::new(&f);
        for line in file.lines().skip(1) {
            match line{
                Ok(l) =>  {
                    let parts = l.split("\t").collect::<Vec<_>>();
                    let sol = GoodSolution{
                        a_nm: parts[0].to_owned(),
                        b_nm: parts[1].to_owned(),
                        or: if parts[2]=="N" {Orientation::Normal} else {Orientation::Reversed},
                        oha: parts[3].parse().unwrap(),
                        ohb: parts[4].parse().unwrap(),
                        ola: parts[5].parse().unwrap(),
                        olb: parts[6].parse().unwrap(),
                        err: parts[7].parse().unwrap(),
                    };
                    result.insert(sol);
                },
                Err(_) => println!(),
            };
        }
        result
    }

    fn panic_if_solutions_missing(solutions : HashSet<GoodSolution>, should_contain : HashSet<GoodSolution>){
        let mut found : HashSet<&GoodSolution> = HashSet::new();
        let mut missing : HashSet<&GoodSolution> = HashSet::new();
        for x in should_contain.iter(){
            if solutions.contains(&x){
                found.insert(x);
            }else{
                missing.insert(x);
            }
        }
        if missing.len() > 0{
            panic!(format!("MISSING {} / {} solutions\n{:#?}. all solutions {:#?}",
                           missing.len(), should_contain.len(), &missing, &solutions));
        }
    }
}