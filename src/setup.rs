use num_cpus;
use structs::run_config::Config;
use modes::{IsMode, Mode};
use modes;
use std::cmp::{min, max};

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
        (@arg mode: -m --mode +takes_value "Uses the filtering scheme mode given options {valimaki, kucherov}. Modes can also be supplied string arguments i.e. 'kucherov_2'. (Default : kucherov_2")

        (@arg format_line: -f --format_line "The first line of the output file will contain a TSV header line.")
        (@arg reversals: -r --reversals "Enables reversals of input strings")
        (@arg inclusions: -i --inclusions "Enables finding of inclusion overlaps (one string within another)")
        (@arg edit_distance: -e --edit_distance "Uses Levenshtein / edit distance instead of Hamming distance")
        (@arg verbose: -v --verbose +multiple "Prints completed steps of the run process")
        (@arg greedy_output: -g --greedy_output "Threads print solutions to output greedily instead of storing them. Limited duplication may arise")
        (@arg print: -p --print "For each solution printed to file, also prints a rough visualization to stdout (mostly for debugging purposes)")
        (@arg no_n: -n --no_n "Omits N symbol from alphabet saving time. Will remove N symbols from input file (with a warning)")
        (@arg track_progress: -t --track_progress "Prints progress bar for completed tasks and ETA to stdout")
    ).get_matches();

    let worker_threads = match matches.value_of("worker_threads") {
        Some(s) => s.parse().unwrap(),
        None => max(1, num_cpus::get()-1),
    };
    let mode : Mode = match matches.value_of("mode") {
        Some(s) => modes::get_mode(s),
        _ => modes::default_mode(),
    };

    let config = Config{
        //required
        input  :            matches.value_of("IN_PATH").unwrap().to_owned(),
        output :            matches.value_of("OUT_PATH").unwrap().to_owned(),
        err_rate :          matches.value_of("ERR_RATE").unwrap().parse().unwrap(),
        thresh :            matches.value_of("THRESH").unwrap().parse().unwrap(),

        //options
        worker_threads :    worker_threads,
        verbosity:          min(matches.occurrences_of("verbose") as u8, 2),

        //opt-in
        reversals :         if matches.occurrences_of("reversals")        >= 1 {true} else {false},
        inclusions :        if matches.occurrences_of("inclusions")       >= 1 {true} else {false},
        edit_distance :     if matches.occurrences_of("edit_distance")    >= 1 {true} else {false},
        greedy_output:      if matches.occurrences_of("greedy_output")    >= 1 {true} else {false},
        print:              if matches.occurrences_of("print")            >= 1 {true} else {false},
        track_progress:     if matches.occurrences_of("track_progress")   >= 1 {true} else {false},
        format_line:        if matches.occurrences_of("format_line")      >= 1 {true} else {false},

        //opt-out
        n_alphabet :        if matches.occurrences_of("no_n")             == 0 {true} else {false},
    };

    assert!(config.thresh > 0, "ERROR! Threshold value must be strictly larger than 0.");
    assert!(config.err_rate >= 0.0 && config.err_rate < 1.0, "ERROR! Error rate limit must be non-negative and smaller than 1. 0 <= e < 1.");
    if !config.reversals{
        println!("WARNING! Reversals are NOT enabled by default. Run with -r flag to enable reversals.");
    }
    (mode, config)
}
