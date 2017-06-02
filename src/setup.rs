use num_cpus;
use structs::run_config::Config;
use modes::{IsMode, Mode};
use modes;

/*
Using Clap, builds a config struct that contains all the user's input
*/
pub fn parse_run_args() -> (Box<IsMode>, Config) {
    let matches = clap_app!(ASPOPsolver =>
        (version: "1.0")
        (author: "Christopher Esterhuyse <christopher.esterhuyse@gmail.com>")
        (about: "Finds approximate suffix prefix overlaps from a given fasta file")

        (@arg IN_PATH: +required +takes_value "Path to the input fasta file")
        (@arg OUT_PATH: +required +takes_value "Path of desired output file")
        (@arg ERR_RATE: +required +takes_value "The max rate of errors in an overlap")
        (@arg THRESH: +required +takes_value "Shortest allowed length of an overlap")

        (@arg worker_threads: -w --worker_threads +takes_value "Number of worker threads used. Defaults to number of logical cpu cores")
        (@arg mode: -m --mode +takes_value "Uses the filtering scheme mode given options {valimaki2, kucherov[1,2,3,4...]} (Default : kucherov2")

        (@arg reversals: -r --reversals "Enables reversals of input strings")
        (@arg inclusions: -i --inclusions "Enables finding of inclusion overlaps (one string within another)")
        (@arg edit_distance: -e --edit_distance "Uses Levenshtein / edit distance instead of Hamming distance")
        (@arg verbose: -v --verbose "Prints completed steps of the run process")
        (@arg greedy_output: -g --greedy_output "Threads print solutions to output greedily instead of storing them. Limited duplication may arise")
        (@arg print: -p --print "For each solution printed to file, also prints a rough visualization to stdout (mostly for debugging purposes)")
        (@arg no_n: -n --no_n "Omits N symbol from alphabet saving time. Will remove N symbols from input file (with a warning)")
        (@arg track_progress: -t --track_progress "Prints progress bar for completed tasks and ETA to stdout")
    ).get_matches();

    let worker_threads = match matches.value_of("worker_threads") {
        Some(s) => s.parse().unwrap(),
        None => num_cpus::get(),
    };
    let mode : Mode = match matches.value_of("mode") {
        Some(s) => get_mode(s),
        _ => default_mode(),
    };

    let config = Config{
        //required
        input  :            matches.value_of("IN_PATH").unwrap().to_owned(),
        output :            matches.value_of("OUT_PATH").unwrap().to_owned(),
        err_rate :          matches.value_of("ERR_RATE").unwrap().parse().unwrap(),
        thresh :            matches.value_of("THRESH").unwrap().parse().unwrap(),

        //options
        worker_threads :    worker_threads,
        verbosity:          matches.occurrences_of("verbose") as u8,

        //opt-in
        reversals :         if matches.occurrences_of("reversals")        >= 1 {true} else {false},
        inclusions :        if matches.occurrences_of("inclusions")       >= 1 {true} else {false},
        edit_distance :     if matches.occurrences_of("edit_distance")    >= 1 {true} else {false},
        greedy_output:      if matches.occurrences_of("greedy_output")    >= 1 {true} else {false},
        print:              if matches.occurrences_of("print")            >= 1 {true} else {false},
        track_progress:     if matches.occurrences_of("track_progress")   >= 1 {true} else {false},

        //opt-out
        n_alphabet :        if matches.occurrences_of("no_n")             == 0 {true} else {false},
    };
    (mode, config)
}

pub fn default_mode() -> Mode {
    Box::new(modes::kucherov::KucherovMode::new(2))
}

fn get_mode(arg : &str) -> Mode {
    let tokens : Vec<&str> = arg.split('_').collect();
    if tokens.len() > 2 {
        panic!("Too many mode tokens! Tokens are delimited by '_'. See -h for help.");
    } else if tokens.len() == 2 {
        let param : i32 = tokens[1].parse().expect("Couldn't interpret the mode argument as a number. See -h for help.");
        if tokens[0] == "kucherov"{
            assert!(param >= 1);
            Box::new(modes::kucherov::KucherovMode::new(param))
        }else{
            panic!(format!("'{}' couldn't be interpreted as a mode. See -h for help.", arg))
        }
    } else if tokens.len() == 1 {
        if tokens[0] == "valimaki2"{
            Box::new(modes::valimaki2::Valimaki2Mode::new())
        }else{
            panic!(format!("'{}' couldn't be interpreted as a mode. See -h for help.", arg))
        }
    }else{
        panic!("Mode argument expects a value! See -h for help.");
    }
}