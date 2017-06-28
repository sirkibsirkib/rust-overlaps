# Rust implementation for Kucherov's algorithm of an ASPOP solver

## Rust and Cargo

To get started, you need to ensure rust and cargo are installed on your system. They can be found at https://www.rust-lang.org/en-US/.
Next, you simply clone this repo, and open the resulting `/rust_overlaps/` directory in your terminal. All of the source code for the project located in `./src/` as you would expect.

You can use Cargo to compile the source code for a project into a binary native to your system. Use the `--release` argument to ensure the code is compiled _with_ all compiler optimizations if you intend to use the system for realistic applications (It makes a massive difference!)
```
cargo build --release
```

This may take a minute, as it has to fetch and compile all of the dependencies. Once it's finished, there will be a new `./target/release/` directory containing the binary you are after: `rust_overlaps` (or `rust_overlaps.exe` on Windows).

## Using the solver

To get acquainted with the use of the solver, simply run the compiled `rust_overlaps` binary in your terminal with just the argument `--help` (or `-h` for short). This will print some useful information to your console, showing you just about everything you need to know about the arguments you will need and the flags you may want.

For example, this is a typical-looking command.

```
rust_overlaps ./data/viral_data.fasta ./outputs/viral_overlap_solutions.tsv 0.012 80 -r -vv -t -w=10
```
Lets talk about these arguments, as for 90% of executions, just these above will give you everything you need to know
* `./data/viral_data.fasta` is the input path, which expects file in FASTA format.
* `./outputs/viral_overlap_solutions.tsv` is the output path. This file will be created and written by the solver.
* `0.012` this is the _error rate limit_ parameter. Overlap solutions with overlaps containing no more than 0.012 errors per overlapping symbol will be in the output solution set.
* `80` This is the _overlap threshold length_ parameter. No overlaps with both of the two overlap lengths shorter than 80 will be in the output solution set.
* `-r` set an optional flag of finding _reversal_ solutions in addition to normal ones. Eg: input string 'AAAACG' will correspond with reversal 'CGTTTT'.
* `-vv` this use of the verbose flag `-v` twice sets verbosity of the program to 2, (the maximum) which is 0 by default.
* `-t` this flag enables task progress, and will print to terminal how many of the total tasks are completed and an ETA.
* `-w=10` this `-w` flag expects a numeric argument for the desired number of _worker threads_ for the execution, which defaults to `max(1, number_of_logical_cores()-1)` if not specified.

## Output format
The output file will be formatted as a TSV, with one line for the header, which looks like this:
```
idA	idB	O	OHA	OHB	OLA	OLB	K
```
Each subsequent line will contain a single overlap solution. If the program is not run with the _greedy output_ flag `-g`, the solutions will be unique and sorted in a lexicographic ordering of the columns from left to right. For sorting, the values in columns `idA` and `idB` are considered _strings_.
* `idA` The ID of the first 'A' string involved in the overlap as represented in the input FASTA file.
* `idB` Same as above, but for the other 'B' string.
* `O` Orientation of the overlap. This field will always be `N` for 'normal' if the solver is executed without flag `-r`; Otherwise, an 'I' for 'inverted' in this field suggests that the string corresponding to `idB` was reversed in this overlap.
* `OHA` Overhang of the A string on the left. The size of the _prefix_ of the A string not involved in the overlap if positive, and _suffix_ otherwise.
* `OHA` Overhang of the B string on the right. The size of the _suffix_ of the B string not involved in the overlap if positive, and _prefix_ otherwise.
* `OLA` Overlap of A; The length of the substring of A involved in the overlap.
* `OLA` Overlap of B; The length of the substring of B involved in the overlap.
* `K` The _error distance_ between strings A and B. If flag `-e` is used, this is defined as _edit distance_ and _Hamming distance_ otherwise.

The output solutions always guarantee that for each solution, `idA` < `idB` when the IDs are ordered as _strings_. Also, the A string of an overlap is never of _reversed_ orientation.