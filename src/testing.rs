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

//    #[test]
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
            greedy_output:        false,
            time:           false,
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

//    #[test]
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
            time:           false,
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

//    #[test]
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
            time:           false,
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
    fn reverse() {
        let config = Config{
            input  :        "./test_input/reverse.fasta".to_owned(),
            output  :       "./test_output/reverse.txt".to_owned(),
            err_rate :      0.2,
            thresh :        5,
            worker_threads: 1,
            reversals :         true,
            inclusions :    false,
            edit_distance : false,
            verbose :       false,
            greedy_output:  false,
            time:           false,
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