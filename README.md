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

To get acquainted with the use of the solver, simply run the compiled `rust_overlaps` binary in your terminal with just the argument `--help` (or `-h` for short). This will print some useful information to your console, showing you just about everything you need to know about the arguments you will need and the arguments you may want.