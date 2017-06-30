# ASPOP Solver (Rust Implementation)
This git repo contains everything you'll need to use our Rust ASPOP solver, which is ready to go.
It comes with the suffix-filter algorithms and associated schemes according to Välimäki's and Kucherov's papers from 2012 and 2014 respectively. See the section below for extending the solver with your own custom schemes.

## Rust and Cargo

To get started, you need to ensure rust and cargo are installed on your system. They can be found at https://www.rust-lang.org/en-US/.
Next, you simply clone this repo, and open the resulting `/rust_overlaps/` directory in your terminal. All of the source code for the project located in `./src/` as you would expect.

You can use Cargo to compile the source code for a project into a binary native to your system. Use the `--release` argument to ensure the code is compiled _with_ all compiler optimizations if you intend to use the system for realistic applications (It makes a massive difference!)
```
cargo build --release
```

This may take a minute, as it has to fetch and compile all of the dependencies. Once it's finished, there will be a new `./target/release/` directory containing the binary you are after: `rust_overlaps` (or `rust_overlaps.exe` on Windows).

## Using the Solver

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

## Output Format
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

## Custom Filtering and Partitioning Schemes
This solver comes with 2 existing schemes, and defaults to that of Kucherov et al (2014).
However, it was also specifically designed so that adding new schemes would be as easy as possbible. To do this, simply follow these steps:
1. Create your own `struct`. I suggest making a new .rs file in `src/modes/` in the fashion of the existing files such as `src/modes/kucherov.rs`. I suggest using this existing mode as a starting point in general.
2. Have your struct use and implement the `IsMode` trait defined in `src/modes/mod.rs`. This requires that your struct implement the following functions:
    * `filter_func` This function defines the behaviour of your 'filtering scheme'. It will be called repeatedly during text index query searches to check whether nodes in the search tree are permitted to introduce symbol errors. As such, this function is called and expected to return an integer, represeting the _cap_ for how many errors a search node with the given properties (defined by the parameters) is permitted to have.
    * `get_block_lengths` This function defines the behaviour of the 'partition scheme'. For a pattern of given length, it expects a sequence of _lengths_. These will be interpreted as the lenths of partition blocks (left to right) to split the pattern string into. As such, the lengths returned here should all sum to arg `patt_len`.
    * `candidate_condition` This function allows you to optionally inhibit candidate generation for a search node conditionally. i.e. if your function is defined simply as `true` then every node of the search tree will generate candidates for all match locations.
    * `get_fewest_suff_blocks` This function defines which queries NOT to initiate. The pattern will only create query searches for pattern-block-sequence suffixes of this length or more. 
    * `get_guaranteed_extra_blocks` This function is only requried for `testing.rs` and the `cargo test` that runs the code within. It is intended to represent how many 0-error blocks your partition scheme gaurantees for valid pattern prefixes. If you have no intention of using the given tests, feel free to define this function as returning a dummy value.
3. Implement some other functions required by IsMode. Namely `std::fmt::Display` and `std::fmt::Debug`. I suggest you just copy and paste from the Kucherov code and make the necessary changes
4. In `src/setup.rs`, go to approximately line 100, with the `YOUR MODES GO HERE ^^^^` comment. Above you will find more detailed instructions. The purpose of this step is to get the solver to use your Mode struct when the program is started with `-m` and an appropriate parameter. Note that your struct can optionally accept user's input delimited by underscores. For example: `-m=kucherov_2` will use the kucherov mode and pass it one parameter, "2" which the struct's constructor will interpret accordingly.
5. Build your edited rust source code as described in the section above, called "Rust and Cargo". 
5. Whenever you use the compiled solver, be sure to pass flag `-m=???` where "???" is whetever you defined it as in step 4 (conceptually, your solver's name). Don't forget the optional arguments if you need them!
