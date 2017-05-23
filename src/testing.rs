
/*
Some simple tests to check if the system and all its various modes behave correctly and find
the right solutions, using some predetermined inputs in th text_input directory.
*/
#[cfg(test)]
mod tests {
    use super::super::*;
    use std::io::{BufReader, BufRead};

    #[derive (Eq, PartialEq, Hash, Debug)]
    struct GoodSolution{
        a_name : String,
        b_name : String,
        orientation : Orientation,
        overhang_left_a : i32,
        overhang_right_b : i32,
        overlap_a : usize,
        overlap_b : usize,
        errors : u32,
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
            sorted :        false,
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
            sorted :        false,
            time:           false,
            print:          false,
            n_alphabet:     false,
        };
        let maps = prepare::read_and_prepare(&config.input, &config).expect("Couldn't interpret data.");
        solve(&config, &maps);
        let results = read_output(&config.output);
        must_contain(&results, "x", "y", Orientation::Normal, 5, 10, 5, 5, 0);
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
            sorted :        false,
            time:           false,
            print:          false,
            n_alphabet:     false,
        };
        let maps = prepare::read_and_prepare(&config.input, &config).expect("Couldn't interpret data.");
        solve(&config, &maps);
        let results = read_output(&config.output);
        must_contain(&results, "x", "y", Orientation::Normal, 7, 7, 6, 5, 1);
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
                        a_name : parts[0].to_owned(),
                        b_name : parts[1].to_owned(),
                        orientation : if parts[2]=="N" {Orientation::Normal} else {Orientation::Reversed},
                        overhang_left_a : parts[3].parse().unwrap(),
                        overhang_right_b : parts[4].parse().unwrap(),
                        overlap_a : parts[5].parse().unwrap(),
                        overlap_b : parts[6].parse().unwrap(),
                        errors : parts[7].parse().unwrap(),
                    };
                    result.insert(sol);
                },
                Err(_) => println!(),
            };
        }
        result
    }

    fn must_contain(solutions : &HashSet<GoodSolution>,
                    a_name : &str, b_name : &str, orientation : Orientation,
                    overhang_left_a : i32, overhang_right_b : i32,
                    overlap_a : usize, overlap_b : usize, errors : u32) {
        let sol = GoodSolution{
            a_name : a_name.to_owned(),
            b_name : b_name.to_owned(),
            orientation : orientation,
            overhang_left_a : overhang_left_a,
            overhang_right_b : overhang_right_b,
            overlap_a : overlap_a,
            overlap_b : overlap_b,
            errors : errors,
        };
        if !solutions.contains(&sol){
            panic!(format!("Solution missing! {:#?}\nSolutions : {:#?}", &sol, &solutions));
        }
    }
}