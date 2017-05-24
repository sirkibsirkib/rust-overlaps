
use structs::run_config::Config;

/*
Using Clap, builds a config struct that contains all the user's input
*/
pub fn parse_run_args() -> Config{
    let matches = clap_app!(ASPOPsolver =>
        (version: "1.0")
        (author: "Christopher Esterhuyse <christopher.esterhuyse@gmail.com>")
        (about: "Finds approximate suffix prefix overlaps from a given fasta file")
        (@arg IN_PATH: +required +takes_value "Path to the input fasta file")
        (@arg OUT_PATH: +required +takes_value "Path of desired output file")
        (@arg ERR_RATE: +required +takes_value "The max rate of errors in an overlap")
        (@arg THRESH: +required +takes_value "Shortest allowed length of an overlap")
        (@arg WORKER_THREADS: +required +takes_value "Number of worker threads used")
        (@arg reversals: -r --reversals "Enables reversals of input strings")
        (@arg inclusions: -i --inclusions "Enables finding of inclusion overlaps (one string within another)")
        (@arg edit_distance: -e --edit_distance "Uses Levenshtein / edit distance instead of Hamming /  distance")
        (@arg verbose: -v --verbose "Prints completed steps of the run process")
        (@arg greedy_output: -g --greedy_output "Threads print solutions to output greedily. Saves time and space but outputs will not be sorted and MAY have the same solutions x2 in some cases")
        (@arg time: -t --time "Records time taken to finish working.")
        (@arg print: -p --print "For each solution printed to file, also prints a rough visualization to stdout (mostly for debugging purposes).")
        (@arg no_n: -n --no_n "Omits 'N' from the alphabet and cleans it from input strings, increasing run speed. Will not consider N when searching.")
    ).get_matches();

    Config{
        input  : matches.value_of("IN_PATH").unwrap().to_owned(),
        output : matches.value_of("OUT_PATH").unwrap().to_owned(),
        err_rate : matches.value_of("ERR_RATE").unwrap().parse().unwrap(),
        thresh : matches.value_of("THRESH").unwrap().parse().unwrap(),
        worker_threads: matches.value_of("WORKER_THREADS").unwrap().parse().unwrap(),
        reversals :     if matches.occurrences_of("reversals")     >= 1 {true} else {false},
        inclusions :    if matches.occurrences_of("inclusions")    >= 1 {true} else {false},
        edit_distance : if matches.occurrences_of("edit_distance") >= 1 {true} else {false},
        verbose :       if matches.occurrences_of("verbose")       >= 1 {true} else {false},
        greedy_output:  if matches.occurrences_of("greedy_output") >= 1 {true} else {false},
        time:           if matches.occurrences_of("time")          >= 1 {true} else {false},
        print:          if matches.occurrences_of("print")         >= 1 {true} else {false},
        n_alphabet :    if matches.occurrences_of("no_n")          == 0 {true} else {false},
    }
}